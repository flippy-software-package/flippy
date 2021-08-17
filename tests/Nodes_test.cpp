#include "external/catch.hpp"
#include "json.hpp"
#include <iostream>
#include "utilities/debug_utils.h"
#include "vec3.hpp"
#include "Nodes.h"

using json = nlohmann::json;
using namespace fp;
const std::vector<std::vector<int>> ICOSA_NN_IDS = {
	{4,3,2,1,5},
	{7,6,2,5,0},
	{8,7,3,1,0},
	{9,8,4,2,0},
	{9,10,3,5,0},
	{6,10,4,1,0},
	{11,7,10,1,5},
	{11,8,6,2,1},
	{11,9,7,3,2},
	{11,8,10,4,3},
	{11,9,6,4,5},
	{9,8,7,6,10}
};
const std::vector<std::vector<double>> ICOSA_POSES = {
	{0.0,0.0,100.0},
	{89.44271909999158,0.0,44.721359549995796},
	{27.639320225002102,85.06508083520399,44.721359549995796},
	{-72.36067977499789,52.57311121191337,44.7213595499958},
	{-72.3606797749979,-52.57311121191336,44.7213595499958},
	{27.639320225002088,-85.065080835204,44.7213595499958},
	{72.36067977499789,-52.57311121191336,-44.72135954999579},
	{72.36067977499789,52.57311121191336,-44.72135954999579},
	{-27.639320225002095,85.06508083520399,-44.72135954999579},
	{-89.44271909999158,1.0953573965284053e-14,-44.72135954999579},
	{-27.639320225002113,-85.06508083520399,-44.72135954999579},
	{1.2246467991473532e-14,0.0,-100.0}
};

json const ICOSA_DATA =
	R"({
	  "0":	{"nn_ids": [4,3,2,1,5],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.0,0.0,100.0]},
	  "1":  {"nn_ids": [7,6,2,5,0],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [89.44271909999158,0.0,44.721359549995796]},
	  "2":  {"nn_ids": [8,7,3,1,0],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [27.639320225002102,85.06508083520399,44.721359549995796]},
	  "3":  {"nn_ids": [9,8,4,2,0],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-72.36067977499789,52.57311121191337,44.7213595499958]},
	  "4":  {"nn_ids": [9,10,3,5,0],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-72.3606797749979,-52.57311121191336,44.7213595499958]},
	  "5":  {"nn_ids": [6,10,4,1,0],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [27.639320225002088,-85.065080835204,44.7213595499958]},
	  "6":  {"nn_ids": [11,7,10,1,5], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [72.36067977499789,-52.57311121191336,-44.72135954999579]},
	  "7":  {"nn_ids": [11,8,6,2,1],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [72.36067977499789,52.57311121191336,-44.72135954999579]},
	  "8":  {"nn_ids": [11,9,7,3,2],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-27.639320225002095,85.06508083520399,-44.72135954999579]},
	  "9":  {"nn_ids": [11,8,10,4,3], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-89.44271909999158,1.0953573965284053e-14,-44.72135954999579]},
	  "10": {"nn_ids": [11,9,6,4,5],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-27.639320225002113,-85.06508083520399,-44.72135954999579]},
	  "11": {"nn_ids": [9,8,7,6,10],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [1.2246467991473532e-14,0.0,-100.0]}
  })"_json;


TEST_CASE("Correct Read of Icosa Data"){

  Nodes<double, short> icosa_nodes(ICOSA_DATA, 0);
  for (int i = 0; i < 12; ++i) {
    CHECK(icosa_nodes.data[i].id==i );

    for (int j = 0; j < 3; ++j) {
	  CHECK(icosa_nodes.data[i].pos[j]==ICOSA_DATA[std::to_string(i)]["pos"][j]);
	}
	for (int j = 0; j < 5; ++j) {
	  CHECK(icosa_nodes.data[i].nn_ids[j]==ICOSA_DATA[std::to_string(i)]["nn_ids"][j]);
	}
  }
}

TEST_CASE("pop emplace test"){
    using real = double;
    using idx = short;
    Node<real, idx> single_node{.id=1, .pos={1,1,1}, .nn_ids={1,3,2}, .nn_distances{vec3<real>{1,0,0}, vec3<real>{0,1,0}, vec3<real>{0,0, 1}}};
    SECTION("simple pop test 1"){
        single_node.pop_nn(3);
        auto ids_res = std::vector<idx>{1,2};
        auto dst_res = std::vector<vec3<real>>{vec3<real>{1,0,0}, vec3<real>{0,0, 1}};
        CHECK(single_node.nn_ids==ids_res);
        CHECK(single_node.nn_distances==dst_res);
    }

    SECTION("simple empl test 1"){
        single_node.emplace_nn_id(7,vec3<real>{0,1,1},2);
        auto ids_res = std::vector<idx>{1,3,7, 2};
        auto dst_res = std::vector<vec3<real>>{vec3<real>{1,0,0}, vec3<real>{0,1,0}, vec3<real>{-1,0,0}, vec3<real>{0,0,1} };
        CHECK(single_node.nn_ids==ids_res);
        CHECK(single_node.nn_distances==dst_res);
    }

}


TEST_CASE("get_distance_to test"){
    using real = double;
    using idx = short;
    Node<real, idx> single_node{.id=1, .pos={1,1,1}, .nn_ids={1,3,2}, .nn_distances{vec3<real>{1,0,0}, vec3<real>{0,1,0}, vec3<real>{0,0, 1}}};
    SECTION("simple get test 1"){
        auto exp_dist = vec3<real>{1,0,0};
        CHECK(single_node.get_distance_vector_to(1)==exp_dist);
        exp_dist = vec3<real>{0,1,0};
        CHECK(single_node.get_distance_vector_to(3)==exp_dist);
        exp_dist = vec3<real>{0,0,1};
        CHECK(single_node.get_distance_vector_to(2)==exp_dist);
    }

    SECTION("simple get after empl test 1"){
        single_node.emplace_nn_id(7,vec3<real>{0,1,1},2);
        auto exp_dist = vec3<real>{0,0,1};
        CHECK(single_node.get_distance_vector_to(2)==exp_dist);
    }

    SECTION("simple get after pop test 1"){
        single_node.pop_nn(3);
        auto exp_dist = vec3<real>{0,0,1};
        CHECK(single_node.get_distance_vector_to(2)==exp_dist);
    }
}
