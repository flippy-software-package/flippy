#include "external/catch.hpp"
#include <array>
#include <iostream>
#include "json.hpp"
#include "utilities/debug_utils.h"

#define private public // please don't do this at home!
#include "vec3.hpp"
#include "Nodes.h"

#include "Triangulation.h"

using json = nlohmann::json;
using namespace fp;
json const ICOSA_DATA =
        R"({
	  "0":	{"nn_ids": [4,2,3,1,5],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.0,0.0,100.0]},
	  "1":  {"nn_ids": [7,6,2,5,0],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [89.44271909999158,0.0,44.721359549995796]},
	  "2":  {"nn_ids": [7,8,3,1,0],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [27.639320225002102,85.06508083520399,44.721359549995796]},
	  "3":  {"nn_ids": [9,8,4,0,2],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-72.36067977499789,52.57311121191337,44.7213595499958]},
	  "4":  {"nn_ids": [9,10,3,5,0],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-72.3606797749979,-52.57311121191336,44.7213595499958]},
	  "5":  {"nn_ids": [6,10,4,1,0],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [27.639320225002088,-85.065080835204,44.7213595499958]},
	  "6":  {"nn_ids": [11,7,10,1,5], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [72.36067977499789,-52.57311121191336,-44.72135954999579]},
	  "7":  {"nn_ids": [11,8,6,2,1],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [72.36067977499789,52.57311121191336,-44.72135954999579]},
	  "8":  {"nn_ids": [11,9,7,3,2],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-27.639320225002095,85.06508083520399,-44.72135954999579]},
	  "9":  {"nn_ids": [11,8,10,4,3], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-89.44271909999158,1.0953573965284053e-14,-44.72135954999579]},
	  "10": {"nn_ids": [11,9,6,4,5],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-27.639320225002113,-85.06508083520399,-44.72135954999579]},
	  "11": {"nn_ids": [9,8,7,6,10],  "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [1.2246467991473532e-14,0.0,-100.0]}
  })"_json;

json const DOUBLE_TETRA_DATA =
        R"({
	  "0":	{"nn_ids": [1,2,3],     "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.0,0.0,1.0]},
	  "1":  {"nn_ids": [0,2,3,4,5], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.0,0.0,0.0]},
	  "2":  {"nn_ids": [0,1,3,4],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [1.0,-2.0,0.0]},
	  "3":  {"nn_ids": [0,1,2,4,5], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [1.0,1,0]},
	  "4":  {"nn_ids": [1,2,3,5],   "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-72.3606797749979,-52.57311121191336,44.7213595499958]},
	  "5":  {"nn_ids": [1,3,4],     "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [27.639320225002088,-85.065080835204,44.7213595499958]},
	  "6":  {"nn_ids": [5],         "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [72.36067977499789,-52.57311121191336,-44.72135954999579]}
  })"_json;

// This topology was created by hand. The nn_ids order matters and makes the surface well oriented.
json const STAR_DATA =
        R"({
	  "2":	 {"nn_ids": [4, 3, 6, 7, 9, 10],       "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.0,0.0,0.0]},
	  "3":  {"nn_ids": [8, 6, 2, 4, 5, 1],         "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [1.0,0.0,0.0]},
	  "6":  {"nn_ids": [3, 8, 0, 7, 2],            "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.5, 0.866025, 0.0]},
	  "7":  {"nn_ids": [2, 6, 0, 9],               "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-0.5, 0.866025, 0]},
	  "9": {"nn_ids": [2, 7, 0, 10],               "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-1,0,0]},
	  "10": {"nn_ids": [4, 2, 9, 0],               "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [-0.5, -0.866025, 0]},
	  "4":  {"nn_ids": [3, 2, 10, 0, 5],           "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.5, -0.866025, 0]},
	  "5":  {"nn_ids": [3, 4, 0, 1],               "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [1.5, -0.866025,0]},
	  "1":   {"nn_ids": [3, 5, 0, 8],              "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [2,0,0]},
	  "8": {"nn_ids": [3, 1, 0, 6],                "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [1.5, 0.866025, 0]},
	  "0":   {"nn_ids": [6, 8, 1, 5, 4, 10, 9, 7], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "scaled_curvature_energy": 0, "pos": [0.5,0,-1]}
  })"_json;

TEST_CASE("Propper scaling with initial Radius for Triangulation Instantiator")
{
    double r_init = 10.1;
//  json all_data = json_read("../../init_files/egg_it_0.json");
    Triangulation<double, long> triangulation(ICOSA_DATA, r_init, 0);
    auto target = Approx(r_init).margin(0.001);

    for (long i = 0; i<triangulation.size(); ++i) {
        CHECK((triangulation[i].pos - triangulation.mass_center()).norm()==target);
    }
}

TEST_CASE("Propper move")
{
    double r_init = 1;
    Triangulation<double, long> icosa_triangulation(ICOSA_DATA, r_init, 0);
    auto new_x_target = Approx(0).margin(0.01);
    auto new_y_target = Approx(0).margin(0.01);
    auto new_z_target = Approx(2).margin(0.01);
//    auto new_mc_z_target = Approx(1./12.).margin(0.01);
    vec3<double> displ{0, 0, 1};

    SECTION("Just CHECK if nodes are moved to the right place") {
        icosa_triangulation.move_node(0, displ);
        CHECK(icosa_triangulation[0].pos[0]==new_x_target);
        CHECK(icosa_triangulation[0].pos[1]==new_y_target);
        CHECK(icosa_triangulation[0].pos[2]==new_z_target);
    }

    SECTION("Just CHECK if move_node reverses things correctly") {
        vec3<double> displ2{0.1, -30, 1};
        Triangulation<double, long> copy_icosa_triangulation{icosa_triangulation};
        icosa_triangulation.move_node(0, displ2);
//        auto copy_icosa_triangulation = icosa_triangulation;
        icosa_triangulation.move_node(0, -displ2);
        CHECK((copy_icosa_triangulation.mass_center()-icosa_triangulation.mass_center()).norm()==Approx(0).margin(0.001));
        for(int i=0; i<icosa_triangulation.size();++i) {
            CHECK(copy_icosa_triangulation.nodes_[i].area==Approx(icosa_triangulation.nodes_[i].area).margin(0.01));
            CHECK(copy_icosa_triangulation.nodes_[i].scaled_curvature_energy==Approx(icosa_triangulation.nodes_[i].scaled_curvature_energy).margin(0.01));
            CHECK(copy_icosa_triangulation.nodes_[i].id==icosa_triangulation.nodes_[i].id);
            CHECK(copy_icosa_triangulation.nodes_[i].nn_ids==icosa_triangulation.nodes_[i].nn_ids);
            for(std::size_t j=0; j < copy_icosa_triangulation.nodes_[i].nn_ids.size();++j) {
                CHECK((copy_icosa_triangulation.nodes_[i].nn_distances[j] - icosa_triangulation.nodes_[i].nn_distances[j])
                        .norm()==Approx(0).margin(0.01));
            }
        }
    }

//    SECTION("Just CHECK if unmove_move reverses things correctly") {
//        vec3<double> displ2{0.1, -30, 1};
//        Triangulation<double, long> copy_icosa_triangulation{icosa_triangulation};
//        icosa_triangulation.move_node(0, displ2);
////        auto copy_icosa_triangulation = icosa_triangulation;
//        icosa_triangulation.unmove_node();
//        CHECK((copy_icosa_triangulation.mass_center()-icosa_triangulation.mass_center()).norm()==Approx(0).margin(0.001));
//        for(int i=0; i<icosa_triangulation.size();++i) {
//            CHECK(copy_icosa_triangulation.nodes_[i].area==Approx(icosa_triangulation.nodes_[i].area).margin(0.01));
//            CHECK(copy_icosa_triangulation.nodes_[i].scaled_curvature_energy==Approx(icosa_triangulation.nodes_[i].scaled_curvature_energy).margin(0.01));
//            CHECK(copy_icosa_triangulation.nodes_[i].id==icosa_triangulation.nodes_[i].id);
//            CHECK(copy_icosa_triangulation.nodes_[i].nn_ids==icosa_triangulation.nodes_[i].nn_ids);
//            for(std::size_t j=0; j< copy_icosa_triangulation.nodes_[i].nn_ids.size();++j) {
//                CHECK((copy_icosa_triangulation.nodes_[i].nn_distances[j] - icosa_triangulation.nodes_[i].nn_distances[j])
//                        .norm()==Approx(0).margin(0.01));
//            }
//        }
//    }
}
TEST_CASE("Propper topology change")
{

    SECTION("CHECK two common neighbours and normal common neighbours on icosa examples") {
        double r_init = 1;
        Triangulation<double, long> icosa_triangulation(ICOSA_DATA, r_init);

        std::array<long, 2> two_cnns = icosa_triangulation.two_common_neighbours(11, 1);
        std::vector<long> cnns = icosa_triangulation.common_neighbours(11, 1);
        std::sort(two_cnns.begin(), two_cnns.end());
        std::sort(cnns.begin(), cnns.end());
        CHECK(two_cnns==std::array<long, 2>{6, 7});
        CHECK(cnns==std::vector<long>{6, 7});

        two_cnns = icosa_triangulation.two_common_neighbours(5, 7);
        cnns = icosa_triangulation.common_neighbours(5, 7);
        std::sort(two_cnns.begin(), two_cnns.end());
        std::sort(cnns.begin(), cnns.end());
        CHECK(two_cnns==std::array<long, 2>{1, 6});
        CHECK(cnns==std::vector<long>{1, 6});
    }

    SECTION("CHECK new two common neighbours on STAR_DATA examples") {
        Triangulation<double, long> star_triangulation(STAR_DATA, 0);
        Neighbors<long> cnns = star_triangulation.previous_and_next_neighbour_global_ids(2l, 3l);
        Neighbors<long> cnns_other_way_around = star_triangulation.previous_and_next_neighbour_global_ids(3l, 2l);

        CHECK(cnns.j_m_1==(long long) 4);
        CHECK(cnns.j_p_1==(long long) 6);
        CHECK(cnns_other_way_around.j_m_1==(long long) 6);
        CHECK(cnns_other_way_around.j_p_1==(long long) 4);
    }

//  Triangulation<double, long> tetra_triangulation(DOUBLE_TETRA_DATA);
//  SECTION("CHECK common neighbors for node pairs with more common neighbors") {
//	std::vector<long> cnns = tetra_triangulation.common_neighbours(0, 4);
//	std::sort(cnns.begin(), cnns.end());
//	CHECK(cnns==std::vector<long>{1, 2, 3});
//
//	cnns = tetra_triangulation.common_neighbours(0, 5);
//	std::sort(cnns.begin(), cnns.end());
//	CHECK(cnns==std::vector<long>{1, 3});
//
//	cnns = tetra_triangulation.common_neighbours(0, 6);
//	CHECK(cnns.size()==0);
//  }
}
TEST_CASE("unittest private member functions")
{
    double r_init = 1;
    Triangulation<double, long> icosa_triangulation(ICOSA_DATA, r_init);
//  Triangulation<double, long> tetra_triangulation(DOUBLE_TETRA_DATA);
//  SECTION("unittest all_nn_distance_vectors") {
//	std::vector<vec3<double>> all_nn_d = tetra_triangulation.all_nn_distance_vectors(0l, std::vector{1l, 2l, 3l});
//	vec3<double> expected01 = vec3<double>{0, 0, -1};
//	vec3<double> expected02 = vec3<double>{1, -2, -1};
//	vec3<double> expected03 = vec3<double>{1, 1, -1};
//	CHECK(all_nn_d[0]==expected01);
//	CHECK(all_nn_d[1]==expected02);
//	CHECK(all_nn_d[2]==expected03);
//	auto dist01 = Approx(1).margin(0.001);
//	auto dist02 = Approx(sqrt(6)).margin(0.001);
//	auto dist03 = Approx(sqrt(3)).margin(0.001);
//	CHECK(all_nn_d[0].norm()==dist01);
//	CHECK(all_nn_d[1].norm()==dist02);
//	CHECK(all_nn_d[2].norm()==dist03);
//  }

    SECTION("cos_bond_opposite_angles and cot_alphas_sum") {
        std::vector<double> coses(2);
//	std::tie(coses[0], coses[1]) = tetra_triangulation.cos_bond_opposite_angles(0, 1);
//	std::sort(coses.begin(), coses.end());
//	CHECK(coses[0]==Approx(0.816497).margin(0.0001));
//	CHECK(coses[1]==Approx(0.912871).margin(0.0001));
//	CHECK(tetra_triangulation.cot_alphas_sum(0, 1)==Approx(3.65028).margin(0.0001));
    }

    SECTION("order_nn_ids: simple test on icosa data") {
        Triangulation<double, long long> icosa_triangulation_load(ICOSA_DATA, 0);
        CHECK(icosa_triangulation_load.order_nn_ids(0)==std::vector<long long>{3, 4, 5, 1, 2});
    }

    SECTION("orient_surface_of_a_sphere: simple test on icosa data") {
        Triangulation<double, long long> icosa_triangulation_load(ICOSA_DATA, 0);
        icosa_triangulation_load.orient_surface_of_a_sphere();

        CHECK(icosa_triangulation_load[0].nn_ids==std::vector<long long>{3, 4, 5, 1, 2});
        CHECK(icosa_triangulation_load[1].nn_ids==std::vector<long long>{6, 7, 2, 0, 5});
        CHECK(icosa_triangulation_load[3].nn_ids==std::vector<long long>{8, 9, 4, 0, 2});
        CHECK(icosa_triangulation_load[10].nn_ids==std::vector<long long>{9, 11, 6, 5, 4});

    }
}

TEST_CASE("unittest private static functions")
{
//    SECTION("unittest cos_between_vectors for a few vectors") {
//        auto li1 = vec3<double>{1, 0, 0};
//        auto li2 = vec3<double>{0, -1, 0};
//        double cos_angle = Triangulation<double, int>::cos_between_vectors(li1, li2);
//        auto cos_angle_target = Approx(0).margin(0.001);
//        CHECK(cos_angle==cos_angle_target);
//
//        li1 = vec3<double>{1, 0, 0};
//        li2 = vec3<double>{1, 0, 0};
//        cos_angle = Triangulation<double, int>::cos_between_vectors(li1, li2);
//        cos_angle_target = Approx(1).margin(0.001);
//        CHECK(cos_angle==cos_angle_target);
//
//        li1 = vec3<double>{1, 0, 0};
//        li2 = vec3<double>{-1, 0, 0};
//        cos_angle = Triangulation<double, int>::cos_between_vectors(li1, li2);
//        cos_angle_target = Approx(-1).margin(0.001);
//        CHECK(cos_angle==cos_angle_target);
//
//        li1 = vec3<double>{sqrt(3), 0, 0};
//        li2 = vec3<double>{sqrt(3), 0, 1};
//        double cos_angle30 = Triangulation<double, int>::cos_between_vectors(li1, li2);
//        cos_angle_target = Approx(std::cos(M_PI/6)).margin(0.001);
//        CHECK(cos_angle30==cos_angle_target);
//
//        li1 = vec3<double>{1, sqrt(3), 0};
//        li2 = vec3<double>{1, 0, 0};
//        double cos_angle60 = Triangulation<double, int>::cos_between_vectors(li1, li2);
//        cos_angle_target = Approx(std::cos(M_PI/3)).margin(0.001);
//        CHECK(cos_angle60==cos_angle_target);
//
//        li1 = vec3<double>{0, 1, 0};
//        li2 = vec3<double>{0, 1, 1};
//        double cos_angle45 = Triangulation<double, int>::cos_between_vectors(li1, li2);
//        cos_angle_target = Approx(std::cos(M_PI/4)).margin(0.001);
//        CHECK(cos_angle45==cos_angle_target);
//
//        li1 = vec3<double>{1, 0, 0};
//        li2 = vec3<double>{-1, 1, 0};
//        double cos_angleMin45 = Triangulation<double, int>::cos_between_vectors(li1, li2);
//        cos_angle_target = Approx(std::cos(3*M_PI/4)).margin(0.001);
//        CHECK(cos_angleMin45==cos_angle_target);
//    }

//    SECTION("unittest mixed_area for a few triangles") {
//        auto li1 = vec3<double>{1, 0, 0};
//        auto li2 = vec3<double>{0, 1, 0};
//        double triangle_area_at_node_orth = li1.cross(li2).norm()/2;
//        double mixed_area = Triangulation<double, int>::mixed_area(li1, li2, triangle_area_at_node_orth);
//        CHECK(mixed_area==Approx(triangle_area_at_node_orth/2).margin(0.001));
//
//        li1 = vec3<double>{1, 0, 0};
//        li2 = vec3<double>{-1, 1, 0};
//        double triangle_area_at_node_obtuse = li1.cross(li2).norm()/2;
//        mixed_area = Triangulation<double, int>::mixed_area(li1, li2, triangle_area_at_node_obtuse);
//        CHECK(mixed_area==Approx(triangle_area_at_node_obtuse/2).margin(0.001));
//
//        li1 = vec3<double>{1, 0, 0};
//        li2 = vec3<double>{2, 1, 0};
//        double triangle_area_not_at_node_obtuse = li1.cross(li2).norm()/2;
//        mixed_area = Triangulation<double, int>::mixed_area(li1, li2, triangle_area_not_at_node_obtuse);
//        CHECK(mixed_area==Approx(triangle_area_not_at_node_obtuse/4).margin(0.001));
//
//        li1 = vec3<double>{1, 0, 0};
//        li2 = vec3<double>{.8, 0.5, 0};
//        double triangle_area_not_obtuse = li1.cross(li2).norm()/2;
//        mixed_area = Triangulation<double, int>::mixed_area(li1, li2, triangle_area_not_obtuse);
//        CHECK(mixed_area==Approx(0.070025070025).margin(0.0000001));
//    }

//    SECTION("Partial voronoi area simple test") {
//        auto li1 = vec3<double>{1, 0, 0};
//        auto li2 = vec3<double>{.8, 0.5, 0};
//        auto outward_normal = li1.cross(li2);
//        auto[mixed_area, face_normal] = Triangulation<double, int>::partial_voronoi_area_and_volume_of_node(li1, li2);
//        CHECK(mixed_area==Approx(0.070025070025).margin(0.0000001));
//        CHECK(outward_normal.dot(face_normal)>0);
//    }

//    SECTION("cos_to_cot simple test") {
//        double cos = -0.666276;
//        CHECK(Triangulation<double, int>::cos_to_cot(cos)==Approx(-0.893484).margin(0.0000001));
//
//        auto random_angle = (double) (M_PI/2);
//        auto std_cos = (double) std::cos(random_angle);
//        auto std_cot = (double) (std_cos/std::sin(random_angle));
//
//        CHECK(Triangulation<double, int>::cos_to_cot(std_cos)==Approx(std_cot).margin(0.0000001));
//
//        random_angle = (double) (M_PI/4);
//        std_cos = (double) std::cos(random_angle);
//        std_cot = (double) (std_cos/std::sin(random_angle));
//
//        CHECK(Triangulation<double, int>::cos_to_cot(std_cos)==Approx(std_cot).margin(0.0000001));
//
//        random_angle = (double) (3*M_PI/4);
//        std_cos = (double) std::cos(random_angle);
//        std_cot = (double) (std_cos/std::sin(random_angle));
//
//        CHECK(Triangulation<double, int>::cos_to_cot(std_cos)==Approx(std_cot).margin(0.0000001));
//
//        random_angle = (double) (rand()%3);
//        std_cos = (double) std::cos(random_angle);
//        std_cot = (double) (std_cos/std::sin(random_angle));
//
//        CHECK(Triangulation<double, int>::cos_to_cot(std_cos)==Approx(std_cot).margin(0.0000001));
//
//    }

}