#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <iostream>
#include "flippy.hpp"
#include "test_helper.h"

using namespace fp;

fp::Json const ICOSA_DATA =
	R"({
	  "0":	{"nn_ids": [4,3,2,1,5],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.0,0.0,100.0]},
	  "1":  {"nn_ids": [7,6,2,5,0],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [89.44271909999158,0.0,44.721359549995796]},
	  "2":  {"nn_ids": [8,7,3,1,0],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [27.639320225002102,85.06508083520399,44.721359549995796]},
	  "3":  {"nn_ids": [9,8,4,2,0],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-72.36067977499789,52.57311121191337,44.7213595499958]},
	  "4":  {"nn_ids": [9,10,3,5,0],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-72.3606797749979,-52.57311121191336,44.7213595499958]},
	  "5":  {"nn_ids": [6,10,4,1,0],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [27.639320225002088,-85.065080835204,44.7213595499958]},
	  "6":  {"nn_ids": [11,7,10,1,5], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [72.36067977499789,-52.57311121191336,-44.72135954999579]},
	  "7":  {"nn_ids": [11,8,6,2,1],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [72.36067977499789,52.57311121191336,-44.72135954999579]},
	  "8":  {"nn_ids": [11,9,7,3,2],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-27.639320225002095,85.06508083520399,-44.72135954999579]},
	  "9":  {"nn_ids": [11,8,10,4,3], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-89.44271909999158,1.0953573965284053e-14,-44.72135954999579]},
	  "10": {"nn_ids": [11,9,6,4,5],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-27.639320225002113,-85.06508083520399,-44.72135954999579]},
	  "11": {"nn_ids": [9,8,7,6,10],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [1.2246467991473532e-14,0.0,-100.0]}
  })"_json;


TEST_CASE("Correct Read of Icosa Data"){

  Nodes icosa_nodes(ICOSA_DATA);
  for (Index i = 0; i < 12; ++i) {
    CHECK(icosa_nodes.data[i].id==i );

    for (Index j = 0; j < 3; ++j) {
	  CHECK_THAT(icosa_nodes.data[i].pos[j], ish(ICOSA_DATA[std::to_string(i)]["pos"][j]));
	}
	for (Index j = 0; j < 5; ++j) {
	  CHECK_THAT(icosa_nodes.data[i].nn_ids[j], ish(ICOSA_DATA[std::to_string(i)]["nn_ids"][j]));
	}
  }
}

TEST_CASE("pop emplace test"){
    Node single_node{.id=1, .area=0, .volume=0, .pos={1,1,1},
                                .curvature_vec=vec3<Real>{0,0, 0}, .nn_ids={1,3,2},
                                .nn_distances{vec3<Real>{1,0,0}, vec3<Real>{0,1,0}, vec3<Real>{0,0, 1}},
                                .verlet_list={}};
    SECTION("simple pop test 1"){
        single_node.pop_nn(3);
        auto ids_res = std::vector<Index>{1,2};
        auto dst_res = std::vector<vec3<Real>>{vec3<Real>{1,0,0}, vec3<Real>{0,0, 1}};
        CHECK(single_node.nn_ids==ids_res);
        CHECK(single_node.nn_distances==dst_res);
    }

    SECTION("simple empl test 1"){
        single_node.emplace_nn_id(7,vec3<Real>{0,1,1},2);
        auto ids_res = std::vector<Index>{1,3,7, 2};
        auto dst_res = std::vector<vec3<Real>>{vec3<Real>{1,0,0}, vec3<Real>{0,1,0}, vec3<Real>{-1,0,0}, vec3<Real>{0,0,1} };
        CHECK(single_node.nn_ids==ids_res);
        CHECK(single_node.nn_distances==dst_res);
    }

}


TEST_CASE("get_distance_to test"){
    Node single_node{.id=1, .area=0, .volume=0, .pos={1,1,1},
                                .curvature_vec=vec3<Real>{0,0, 0}, .nn_ids={1,3,2},
                                .nn_distances{vec3<Real>{1,0,0}, vec3<Real>{0,1,0}, vec3<Real>{0,0, 1}},
                                .verlet_list={}};
    SECTION("simple get test 1"){
        auto exp_dist = vec3<Real>{1,0,0};
        CHECK(single_node.get_distance_vector_to(1)==exp_dist);
        exp_dist = vec3<Real>{0,1,0};
        CHECK(single_node.get_distance_vector_to(3)==exp_dist);
        exp_dist = vec3<Real>{0,0,1};
        CHECK(single_node.get_distance_vector_to(2)==exp_dist);
    }

    SECTION("simple get after empl test 1"){
        single_node.emplace_nn_id(7,vec3<Real>{0,1,1},2);
        auto exp_dist = vec3<Real>{0,0,1};
        CHECK(single_node.get_distance_vector_to(2)==exp_dist);
    }

    SECTION("simple get after pop test 1"){
        single_node.pop_nn(3);
        auto exp_dist = vec3<Real>{0,0,1};
        CHECK(single_node.get_distance_vector_to(2)==exp_dist);
    }
}

TEST_CASE("getter and setter tests for Nodes"){
    Nodes icosa_nodes(ICOSA_DATA);

    SECTION("nn_ids"){
        CHECK(icosa_nodes[0].nn_ids ==std::vector<Index>{4,3,2,1,5});
        CHECK(icosa_nodes[1].nn_ids ==std::vector<Index>{7,6,2,5,0});
        CHECK(icosa_nodes[2].nn_ids ==std::vector<Index>{8,7,3,1,0});
        CHECK(icosa_nodes[3].nn_ids ==std::vector<Index>{9,8,4,2,0});
        CHECK(icosa_nodes[4].nn_ids ==std::vector<Index>{9,10,3,5,0});
        CHECK(icosa_nodes[5].nn_ids ==std::vector<Index>{6,10,4,1,0});
        CHECK(icosa_nodes[6].nn_ids ==std::vector<Index>{11,7,10,1,5});
        CHECK(icosa_nodes[7].nn_ids ==std::vector<Index>{11,8,6,2,1});
        CHECK(icosa_nodes[8].nn_ids ==std::vector<Index>{11,9,7,3,2});
        CHECK(icosa_nodes[9].nn_ids ==std::vector<Index>{11,8,10,4,3});
        CHECK(icosa_nodes[10].nn_ids==std::vector<Index>{11,9,6,4,5});
        CHECK(icosa_nodes[11].nn_ids==std::vector<Index>{9,8,7,6,10});
    }


    SECTION("pos"){
        CHECK(icosa_nodes[0].pos ==vec3<Real>{0.0,0.0,100.0});
        CHECK(icosa_nodes[1].pos ==vec3<Real>{89.44271909999158_r, 0.0, 44.721359549995796_r});
        CHECK(icosa_nodes[2].pos ==vec3<Real>{27.639320225002102_r, 85.06508083520399_r, 44.721359549995796_r});
        CHECK(icosa_nodes[3].pos ==vec3<Real>{-72.36067977499789_r, 52.57311121191337_r, 44.7213595499958_r});
        CHECK(icosa_nodes[4].pos ==vec3<Real>{-72.3606797749979_r, -52.57311121191336_r, 44.7213595499958_r});
        CHECK(icosa_nodes[5].pos ==vec3<Real>{27.639320225002088_r, -85.065080835204_r, 44.7213595499958_r});
        CHECK(icosa_nodes[6].pos ==vec3<Real>{72.36067977499789_r, -52.57311121191336_r, -44.72135954999579_r});
        CHECK(icosa_nodes[7].pos ==vec3<Real>{72.36067977499789_r, 52.57311121191336_r, -44.72135954999579_r});
        CHECK(icosa_nodes[8].pos ==vec3<Real>{-27.639320225002095_r, 85.06508083520399_r, -44.72135954999579_r});
        CHECK(icosa_nodes[9].pos ==vec3<Real>{-89.44271909999158_r, 1.0953573965284053e-14_r, -44.72135954999579_r});
        CHECK(icosa_nodes[10].pos==vec3<Real>{-27.639320225002113_r, -85.06508083520399_r, -44.72135954999579_r});
        CHECK(icosa_nodes[11].pos==vec3<Real>{1.2246467991473532e-14_r, 0.0,-100.0});
    }

}