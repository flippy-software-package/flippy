#include "external/catch.hpp"
#include <array>
#include <iostream>
#include "flippy.hpp"

template<typename Index> std::string edge_namer(Index a, Index b){
    if(a<b){
        return std::to_string(a)+"_"+std::to_string(b);
    }else{
        return std::to_string(b)+"_"+std::to_string(a);
    }
}

template<typename Index> std::string face_namer(Index a, Index b, Index c){
    std::array<Index,3> arr{a,b,c};
    std::sort(arr.begin(), arr.end());
    return std::to_string(arr[0]) + "_"
        + std::to_string(arr[1]) + "_"
        + std::to_string(arr[2]);
}

template<typename Index>  std::array<Index, 2> get_two_common_neighbours(std::vector<Index> nn_arr_0, std::vector<Index> nn_arr_1){
    std::array<Index, 2> res{-1, -1};
    for (auto res_p = res.begin(); auto n0_nn_id: nn_arr_0) {
        if (res_p==res.end()) { break; }
        else {
            if (fp::is_member(nn_arr_1, n0_nn_id)) {
                *res_p = n0_nn_id;
                ++res_p;
            }
        }
    }
    return res;
}

TEST_CASE("correct euler number up to nIter=31 count"){
    std::unordered_set<std::string> face_name_hash;
    std::unordered_set<std::string> edge_name_hash;
    std::string edge_name, face_name_0, face_name_1;
    for(short nIter=0; nIter<=31;++nIter){
        fp::Triangulation<float, short, fp::SPHERICAL_TRIANGULATION> trg(nIter, 1.f, 0.f);
        for (auto const& node: trg.nodes()) {
            for(auto nn_id: node.nn_ids){
                edge_name = edge_namer(node.id, nn_id);
                auto cnns = get_two_common_neighbours(node.nn_ids, trg.nodes().nn_ids(nn_id));
                face_name_0 = face_namer(node.id, nn_id, cnns[0]);
                face_name_1 = face_namer(node.id, nn_id, cnns[1]);
                edge_name_hash.insert(edge_name);
                face_name_hash.insert(face_name_0);
                face_name_hash.insert(face_name_1);
            }
        }
        size_t node_count =  trg.nodes().size();
        size_t edge_count =  edge_name_hash.size();
        size_t face_count =  face_name_hash.size();

        SECTION("Euler characteristic"){
            CHECK(node_count - edge_count + face_count == 2);
        }
        SECTION("face edge relation for triangulations"){
            CHECK(edge_count == 3*face_count/2);
        }
    }

}