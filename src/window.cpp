/*!
 * @file window.cpp
 *
 * @brief Window class source file
 */

#include <algorithm>

#include "window.hpp"

#include "abpoa.h"

namespace racon {

std::shared_ptr<Window> createWindow(uint64_t id, uint32_t rank, WindowType type,
    const char* backbone, uint32_t backbone_length, const char* quality,
    uint32_t quality_length) {

    if (backbone_length == 0 || backbone_length != quality_length) {
        fprintf(stderr, "[racon::createWindow] error: "
            "empty backbone sequence/unequal quality length!\n");
        exit(1);
    }

    return std::shared_ptr<Window>(new Window(id, rank, type, backbone,
        backbone_length, quality, quality_length));
}

Window::Window(uint64_t id, uint32_t rank, WindowType type, const char* backbone,
    uint32_t backbone_length, const char* quality, uint32_t quality_length)
        : id_(id), rank_(rank), type_(type), consensus_(), sequences_(),
        qualities_(), positions_() {

    sequences_.emplace_back(backbone, backbone_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(0, 0);
}

Window::~Window() {
}

void Window::add_layer(const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length, uint32_t begin, uint32_t end) {

    if (sequence_length == 0 || begin == end) {
        return;
    }

    if (quality != nullptr && sequence_length != quality_length) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "unequal quality size!\n");
        exit(1);
    }
    if (begin >= end || begin > sequences_.front().second || end > sequences_.front().second) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "layer begin and end positions are invalid!\n");
        exit(1);
    }

    sequences_.emplace_back(sequence, sequence_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(begin, end);
}

struct AbpoaEngine {
    abpoa_t* ab;
    abpoa_para_t* abpt;
    ~AbpoaEngine() {
        abpoa_free(ab);
        abpoa_free_para(abpt);
    }
};

std::shared_ptr<AbpoaEngine> createAbpoaEngine(int8_t match, int8_t mismatch,
    int8_t gap) {

    // racon passes mismatch/gap as negative scores; abPOA expects positive
    // penalties. gap_open1 = 0 selects linear gap mode (as in spoa kNW).
    auto engine = std::make_shared<AbpoaEngine>();
    engine->ab = abpoa_init();
    engine->abpt = abpoa_init_para();
    engine->abpt->align_mode = ABPOA_GLOBAL_MODE;
    engine->abpt->match = match;
    engine->abpt->mismatch = -mismatch;
    engine->abpt->gap_open1 = 0;
    engine->abpt->gap_ext1 = -gap;
    engine->abpt->gap_open2 = 0;
    engine->abpt->gap_ext2 = -gap;
    engine->abpt->out_cons = 1;
    engine->abpt->out_msa = 0;
    engine->abpt->cons_algrm = ABPOA_HB;
    engine->abpt->sub_aln = 1;
    // experiment knobs (temporary): band width and gap model
    if (getenv("RACON_ABPOA_WB")) engine->abpt->wb = atoi(getenv("RACON_ABPOA_WB"));
    if (getenv("RACON_ABPOA_GAP")) {  // "o1,e1,o2,e2" e.g. affine "4,2,0,0"
        sscanf(getenv("RACON_ABPOA_GAP"), "%d,%d,%d,%d", &engine->abpt->gap_open1,
            &engine->abpt->gap_ext1, &engine->abpt->gap_open2, &engine->abpt->gap_ext2);
    }
    abpoa_post_set_para(engine->abpt);
    return engine;
}

static const uint8_t* nt4_table() {
    static uint8_t table[256];
    static bool initialized = [] {
        for (uint32_t i = 0; i < 256; ++i) {
            table[i] = 4;
        }
        table['A'] = table['a'] = 0;
        table['C'] = table['c'] = 1;
        table['G'] = table['g'] = 2;
        table['T'] = table['t'] = 3;
        return true;
    }();
    (void) initialized;
    return table;
}

bool Window::generate_consensus(std::shared_ptr<AbpoaEngine> alignment_engine,
    bool trim) {

    if (sequences_.size() < 3) {
        consensus_ = std::string(sequences_.front().first, sequences_.front().second);
        return false;
    }

    abpoa_t* ab = alignment_engine->ab;
    abpoa_para_t* abpt = alignment_engine->abpt;
    const uint8_t* nt4 = nt4_table();

    abpoa_reset(ab, abpt, sequences_.front().second);
    ab->abs->n_seq = sequences_.size();

    std::vector<uint8_t> bases;
    std::vector<int32_t> weights;
    auto encode = [&](uint32_t i) -> void {
        bases.resize(sequences_[i].second);
        weights.resize(sequences_[i].second);
        for (uint32_t k = 0; k < sequences_[i].second; ++k) {
            bases[k] = nt4[static_cast<uint8_t>(sequences_[i].first[k])];
            // as in spoa: base weight = phred quality, or 1 without qualities
            weights[k] = qualities_[i].first == nullptr ? 1 :
                static_cast<int32_t>(qualities_[i].first[k]) - 33;
        }
    };

    // backbone goes into an empty graph, so backbone position p is node p + 2
    encode(0);
    abpoa_res_t res;
    res.graph_cigar = nullptr; res.n_cigar = 0;
    abpoa_add_graph_alignment(ab, abpt, bases.data(), weights.data(),
        bases.size(), nullptr, res, 0, sequences_.size(), 1);

    std::vector<uint32_t> rank;
    rank.reserve(sequences_.size());
    for (uint32_t i = 0; i < sequences_.size(); ++i) {
        rank.emplace_back(i);
    }

    std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
        return positions_[lhs].first < positions_[rhs].first; });

    uint32_t offset = 0.01 * sequences_.front().second;
    for (uint32_t j = 1; j < sequences_.size(); ++j) {
        uint32_t i = rank[j];
        encode(i);

        res.graph_cigar = nullptr; res.n_cigar = 0;
        if (positions_[i].first < offset && positions_[i].second >
            sequences_.front().second - offset) {
            abpoa_align_sequence_to_graph(ab, abpt, bases.data(), bases.size(), &res);
            abpoa_add_graph_alignment(ab, abpt, bases.data(), weights.data(),
                bases.size(), nullptr, res, j, sequences_.size(), 1);
        } else {
            int exc_beg, exc_end;
            abpoa_subgraph_nodes(ab, abpt, positions_[i].first + 2,
                positions_[i].second + 2, &exc_beg, &exc_end);
            abpoa_align_sequence_to_subgraph(ab, abpt, exc_beg, exc_end,
                bases.data(), bases.size(), &res);
            abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end, bases.data(),
                weights.data(), bases.size(), nullptr, res, j, sequences_.size(), 0);
        }
        if (res.n_cigar) {
            free(res.graph_cigar);
        }
    }

    abpoa_generate_consensus(ab, abpt);
    abpoa_cons_t* abc = ab->abc;
    if (abc->n_cons == 0 || abc->cons_len[0] == 0) {
        consensus_ = std::string(sequences_.front().first, sequences_.front().second);
        return false;
    }

    static const char kBases[] = "ACGTN";
    consensus_.resize(abc->cons_len[0]);
    std::vector<uint32_t> coverages(abc->cons_len[0]);
    for (int32_t k = 0; k < abc->cons_len[0]; ++k) {
        consensus_[k] = kBases[abc->cons_base[0][k]];
        coverages[k] = abc->cons_cov[0][k];
    }

    if (type_ == WindowType::kTGS && trim) {
        uint32_t average_coverage = (sequences_.size() - 1) / 2;

        int32_t begin = 0, end = consensus_.size() - 1;
        for (; begin < static_cast<int32_t>(consensus_.size()); ++begin) {
            if (coverages[begin] >= average_coverage) {
                break;
            }
        }
        for (; end >= 0; --end) {
            if (coverages[end] >= average_coverage) {
                break;
            }
        }

        if (begin >= end) {
            fprintf(stderr, "[racon::Window::generate_consensus] warning: "
                "contig %lu might be chimeric in window %u!\n", id_, rank_);
        } else {
            consensus_ = consensus_.substr(begin, end - begin + 1);
        }
    }

    return true;
}

}
