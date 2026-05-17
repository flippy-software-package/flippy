#include "flippy.hpp"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <iostream>

namespace {
using fp::operator""_r;

template <typename Index> std::string edge_namer(Index a, Index b) {
    if (a < b) {
        return std::to_string(a) + "_" + std::to_string(b);
    } else {
        return std::to_string(b) + "_" + std::to_string(a);
    }
}

template <typename Index> std::string face_namer(Index a, Index b, Index c) {
    std::array<Index, 3> arr{a, b, c};
    std::sort(arr.begin(), arr.end());
    return std::to_string(arr[0]) + "_" + std::to_string(arr[1]) + "_" + std::to_string(arr[2]);
}

std::array<fp::Index, 2> get_two_common_neighbours(std::vector<fp::Index> nn_arr_0, std::vector<fp::Index> nn_arr_1) {
    static const fp::Index   vln = std::numeric_limits<fp::Index>::max();
    std::array<fp::Index, 2> res{vln, vln};
    for (auto res_p = res.begin(); auto n0_nn_id : nn_arr_0) {
        if (res_p == res.end()) {
            break;
        } else {
            if (fp::is_member(nn_arr_1, n0_nn_id)) {
                *res_p = n0_nn_id;
                ++res_p;
            }
        }
    }
    return res;
}
} // namespace

TEST_CASE("correct euler number up to nIter=31 count") {
    std::unordered_set<std::string> face_name_hash;
    std::unordered_set<std::string> edge_name_hash;

    std::string edge_name;
    std::string face_name_0;
    std::string face_name_1;
    for (std::uint16_t nIter = 0; nIter <= 31; ++nIter) {
        auto trg = fp::Triangulation::make_spherical_triangulation(nIter, 1_r, 0_r);
        for (const auto &node : trg.nodes()) {
            for (auto nn_id : node.nn_ids) {
                edge_name   = edge_namer(node.id, nn_id);
                auto cnns   = get_two_common_neighbours(node.nn_ids, trg[nn_id].nn_ids);
                face_name_0 = face_namer(node.id, nn_id, cnns[0]);
                face_name_1 = face_namer(node.id, nn_id, cnns[1]);
                edge_name_hash.insert(edge_name);
                face_name_hash.insert(face_name_0);
                face_name_hash.insert(face_name_1);
            }
        }
        std::size_t node_count = trg.nodes().size();
        std::size_t edge_count = edge_name_hash.size();
        std::size_t face_count = face_name_hash.size();

        SECTION("Euler characteristic") { CHECK(node_count - edge_count + face_count == 2); }
        SECTION("face edge relation for triangulations") { CHECK(edge_count == 3 * face_count / 2); }
    }
}
