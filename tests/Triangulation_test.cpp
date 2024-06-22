#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <array>
#include <iostream>

#define TESTING_FLIPPY_TRIANGULATION_ndh6jclc0qnp274b = 1
#include "flippy.hpp"
using namespace fp;

static void rescale_triangulation(Real R, Triangulation<SPHERICAL_TRIANGULATION>& tr)
{
    tr.R_initial=R;
    tr.scale_all_nodes_to_R_init();
    tr.orient_surface_of_a_sphere();
    tr.initiate_distance_vectors(); //Todo if this is done before orient surface we can save time
    tr.make_global_geometry();
    tr.make_verlet_list();
}
static constexpr auto max_float = 3.40282347e+38F;

template<typename Index>
std::vector<Index> rotate_left(std::vector<Index> v){
    std::rotate(v.begin(), v.begin() + 1, v.end());
    return v;
}

fp::Json const ICOSA_DATA =
        R"({
	  "0":	{"nn_ids": [4,2,3,1,5],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.0,0.0,100.0]},
	  "1":  {"nn_ids": [7,6,2,5,0],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [89.44271909999158,0.0,44.721359549995796]},
	  "2":  {"nn_ids": [7,8,3,1,0],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [27.639320225002102,85.06508083520399,44.721359549995796]},
	  "3":  {"nn_ids": [9,8,4,0,2],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-72.36067977499789,52.57311121191337,44.7213595499958]},
	  "4":  {"nn_ids": [9,10,3,5,0],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-72.3606797749979,-52.57311121191336,44.7213595499958]},
	  "5":  {"nn_ids": [6,10,4,1,0],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [27.639320225002088,-85.065080835204,44.7213595499958]},
	  "6":  {"nn_ids": [11,7,10,1,5], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [72.36067977499789,-52.57311121191336,-44.72135954999579]},
	  "7":  {"nn_ids": [11,8,6,2,1],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [72.36067977499789,52.57311121191336,-44.72135954999579]},
	  "8":  {"nn_ids": [11,9,7,3,2],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-27.639320225002095,85.06508083520399,-44.72135954999579]},
	  "9":  {"nn_ids": [11,8,10,4,3], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-89.44271909999158,1.0953573965284053e-14,-44.72135954999579]},
	  "10": {"nn_ids": [11,9,6,4,5],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-27.639320225002113,-85.06508083520399,-44.72135954999579]},
	  "11": {"nn_ids": [9,8,7,6,10],  "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [1.2246467991473532e-14,0.0,-100.0]}
  })"_json;

// This topology was created by hand. The nn_ids order matters and makes the surface well oriented.
fp::Json const STAR_DATA =
        R"({
	  "2":	 {"nn_ids": [4, 3, 6, 7, 9, 10],       "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.0,0.0,0.0]},
	  "3":  {"nn_ids": [8, 6, 2, 4, 5, 1],         "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [1.0,0.0,0.0]},
	  "6":  {"nn_ids": [3, 8, 0, 7, 2],            "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.5, 0.866025, 0.0]},
	  "7":  {"nn_ids": [2, 6, 0, 9],               "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-0.5, 0.866025, 0]},
	  "9": {"nn_ids": [2, 7, 0, 10],               "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-1,0,0]},
	  "10": {"nn_ids": [4, 2, 9, 0],               "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [-0.5, -0.866025, 0]},
	  "4":  {"nn_ids": [3, 2, 10, 0, 5],           "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.5, -0.866025, 0]},
	  "5":  {"nn_ids": [3, 4, 0, 1],               "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [1.5, -0.866025,0]},
	  "1":   {"nn_ids": [3, 5, 0, 8],              "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [2,0,0]},
	  "8": {"nn_ids": [3, 1, 0, 6],                "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [1.5, 0.866025, 0]},
	  "0":   {"nn_ids": [6, 8, 1, 5, 4, 10, 9, 7], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.5,0,-1]}
  })"_json;

void radius_scaling_test(Triangulation<SPHERICAL_TRIANGULATION> const& triangulation, Real r_init){
    auto target = Catch::Approx(r_init).margin(0.001);
    auto mc = triangulation.calculate_mass_center();
    for (Index i = 0; i<triangulation.size(); ++i) {
        CHECK((triangulation[i].pos - mc).norm()==target);
    }

}

TEST_CASE("Proper scaling with initial Radius for Triangulation Instantiator")
{
    typedef Triangulation<SPHERICAL_TRIANGULATION> SphericalTriangulation;

    Real r_init = 10.1_r;
    SECTION("loading data instantiation"){
        SphericalTriangulation triangulation(ICOSA_DATA, 0);
        rescale_triangulation(r_init, triangulation);
        radius_scaling_test(triangulation, r_init);
    }
    SECTION("spherical triangulator data instantiation"){
        SphericalTriangulation triangulation(0, r_init, 0);
        radius_scaling_test(triangulation, r_init);
    }

}

TEST_CASE("Triangulator Instantiation: correct  node count and global geometry"){
    for (unsigned int nIter = 0; nIter<16; ++nIter) {
        unsigned int nBulk = nIter*(nIter-1)/2;
        unsigned int expected_node_count =
                fp::implementation::IcosahedronSubTriangulation<Real, Index>::N_ICOSA_NODEs
                    + fp::implementation::IcosahedronSubTriangulation<Real, Index>::N_ICOSA_EDGEs*nIter
                    + fp::implementation::IcosahedronSubTriangulation<Real, Index>::N_ICOSA_FACEs * nBulk;
        Triangulation<SPHERICAL_TRIANGULATION> trg(nIter, 1., 0.);
        CHECK(trg.size()==expected_node_count);
        if(nIter>7){
            Real niter_inv = 1_r/(static_cast<Real>(nIter));
            Real precision = 0.14_r*niter_inv;
            auto unit_sphere_volume = Catch::Approx(4_r*PI/3_r).epsilon(precision);
            auto unit_sphere_area = Catch::Approx(4_r*PI).epsilon(precision);
//            auto eight_PI = Catch::Approx(8.*PI).epsilon(precision);
            CHECK(trg.global_geometry().volume==unit_sphere_volume);
            CHECK(trg.global_geometry().area==unit_sphere_area);
//            CHECK(trg.global_geometry().unit_bending_energy==eight_PI);
        }

    }
}

TEST_CASE("Proper move")
{
    Real r_init = 1_r;
    using SphericalTriangulation = Triangulation<SPHERICAL_TRIANGULATION>;
    SphericalTriangulation icosa_triangulation (ICOSA_DATA, 0);
    rescale_triangulation(r_init, icosa_triangulation);
    auto new_x_target = Catch::Approx(0).margin(0.01);
    auto new_y_target = Catch::Approx(0).margin(0.01);
    auto new_z_target = Catch::Approx(2).margin(0.01);
    vec3<Real> displ{0, 0, 1};

    SECTION("Just CHECK if nodes are moved to the right place") {
        icosa_triangulation.move_node(0, displ);
        CHECK(icosa_triangulation[0].pos[0]==new_x_target);
        CHECK(icosa_triangulation[0].pos[1]==new_y_target);
        CHECK(icosa_triangulation[0].pos[2]==new_z_target);
    }

    SECTION("Just CHECK if move_node reverses things correctly") {
        vec3<Real> displ2{0.1_r, -30_r, 1_r};
        SphericalTriangulation copy_icosa_triangulation{icosa_triangulation};
        icosa_triangulation.move_node(0, displ2);
        icosa_triangulation.move_node(0, -displ2);
        vec3<Real> copy_mc = copy_icosa_triangulation.calculate_mass_center();
        vec3<Real> mc = icosa_triangulation.calculate_mass_center();
        CHECK((copy_mc-mc).norm()
                ==Catch::Approx(0).margin(0.001));
        for (Index i = 0; i<icosa_triangulation.size(); ++i) {
            CHECK(copy_icosa_triangulation.nodes_[i].area==Catch::Approx(icosa_triangulation.nodes_[i].area).margin(0.01));
//            CHECK(copy_icosa_triangulation.nodes_[i].unit_bending_energy
//                    ==Catch::Approx(icosa_triangulation.nodes_[i].unit_bending_energy).margin(0.01));
            CHECK(copy_icosa_triangulation.nodes_[i].id==icosa_triangulation.nodes_[i].id);
            CHECK(copy_icosa_triangulation.nodes_[i].nn_ids==icosa_triangulation.nodes_[i].nn_ids);
            for (std::size_t j = 0; j<copy_icosa_triangulation.nodes_[i].nn_ids.size(); ++j) {
                CHECK((copy_icosa_triangulation.nodes_[i].nn_distances[j]
                        - icosa_triangulation.nodes_[i].nn_distances[j]).norm()==Catch::Approx(0).margin(0.01));
            }
        }
    }

}

TEST_CASE("Proper topology change")
{

    SECTION("flip_bond unit test"){
        Triangulation<SPHERICAL_TRIANGULATION> icosa(ICOSA_DATA, 0);
        icosa.orient_surface_of_a_sphere();
        auto bfd = icosa.flip_bond(1,2, 0, max_float);
        CHECK(bfd.flipped==true);
        CHECK(icosa[1].nn_ids==std::vector<Index>{6,7,0,5});
        CHECK(icosa[2].nn_ids==std::vector<Index>{3,0,7,8});
        CHECK(icosa[0].nn_ids==std::vector<Index>{3,4,5,1,7,2});
        CHECK(icosa[7].nn_ids==std::vector<Index>{2,0,1,6,11,8});
    }

    SECTION("Property check: unflip reverses flip"){
        Real r_init = 1_r;
        Real small_number = 0.0001_r;
        Triangulation<SPHERICAL_TRIANGULATION> sphere(10, r_init, 0);
        std::random_device rd;
        std::mt19937_64 rng(rd());
        std::uniform_int_distribution<Index> rid(0, sphere.size()-1);
        for(unsigned int repeat=0; repeat<10; ++repeat) {
            Index rand_node_id = rid(rng);
            auto nn_ids = sphere[rand_node_id].nn_ids;
            std::uniform_int_distribution<Index> r_nn_id(0, static_cast<Index>(nn_ids.size())-1);
            Index rand_nn_id = nn_ids[r_nn_id(rng)];
            auto global_geometry_ORIGINAL = sphere.global_geometry();
            auto bfd = sphere.flip_bond(rand_node_id, rand_nn_id, 0, max_float);
            sphere.unflip_bond(rand_node_id, rand_nn_id, bfd);
            auto global_geometry_AFTER_FLIP_BACKFLIP = sphere.global_geometry();
            CHECK(bfd.flipped==true);
            CHECK(Catch::Approx(global_geometry_ORIGINAL.area).margin(small_number)==global_geometry_AFTER_FLIP_BACKFLIP.area);
            CHECK(Catch::Approx(global_geometry_ORIGINAL.volume).margin(small_number)==global_geometry_AFTER_FLIP_BACKFLIP.volume);
            // todo global energy calculation
//            CHECK(Catch::Approx(global_geometry_ORIGINAL.unit_bending_energy).margin(small_number)==global_geometry_AFTER_FLIP_BACKFLIP.unit_bending_energy);
            bool nn_ids_are_directly_equal_or_equal_after_odd_perm = (nn_ids==sphere[rand_node_id].nn_ids)||(nn_ids==rotate_left(sphere[rand_node_id].nn_ids));
            CHECK(nn_ids_are_directly_equal_or_equal_after_odd_perm);
        }
    }
    SECTION("CHECK two common neighbours and normal common neighbours on icosa examples") {
        Triangulation<SPHERICAL_TRIANGULATION> icosa_triangulation(ICOSA_DATA, 0);

        std::array<Index, 2> two_cnns = icosa_triangulation.two_common_neighbours(11, 1);
        std::vector<Index> cnns = icosa_triangulation.common_neighbours(11, 1);
        std::sort(two_cnns.begin(), two_cnns.end());
        std::sort(cnns.begin(), cnns.end());
        CHECK(two_cnns==std::array<Index, 2>{6, 7});
        CHECK(cnns==std::vector<Index>{6, 7});

        two_cnns = icosa_triangulation.two_common_neighbours(5, 7);
        cnns = icosa_triangulation.common_neighbours(5, 7);
        std::sort(two_cnns.begin(), two_cnns.end());
        std::sort(cnns.begin(), cnns.end());
        CHECK(two_cnns==std::array<Index, 2>{1, 6});
        CHECK(cnns==std::vector<Index>{1, 6});
    }

    SECTION("CHECK new two common neighbours on STAR_DATA examples") {
        Triangulation<SPHERICAL_TRIANGULATION> star_triangulation(STAR_DATA, 0);
        Neighbors cnns = star_triangulation.previous_and_next_neighbour_global_ids(2l, 3l);
        Neighbors cnns_other_way_around = star_triangulation.previous_and_next_neighbour_global_ids(3l, 2l);

        CHECK(cnns.j_m_1==4);
        CHECK(cnns.j_p_1==6);
        CHECK(cnns_other_way_around.j_m_1==6);
        CHECK(cnns_other_way_around.j_p_1==4);
    }

    SECTION("property test: common_neighbours and two_common_neighbours should give same results apart from ordering") {
        auto all_data = json_read("../../tests/init_files/egg_init.json");
        Json nodes = all_data["nodes"];
        Triangulation<SPHERICAL_TRIANGULATION> sphere(nodes, 0);
        std::array<Index, 2> m1m1{VERY_LARGE_NUMBER_, VERY_LARGE_NUMBER_};
        for (Index i = 0; i<sphere.nodes().size() - 1; ++i) {
            for (Index j = i + 1; j<sphere.nodes().size(); ++j) {
                std::vector<Index> cnns = sphere.common_neighbours(i, j);
                std::array<Index, 2> cnns2_arr = sphere.two_common_neighbours(i, j);
                std::vector<Index> cnns2{cnns2_arr[0], cnns2_arr[1]};
                std::sort(cnns.begin(), cnns.end());
                std::sort(cnns2.begin(), cnns2.end());
                if (cnns.empty()) {
                    CHECK(cnns2_arr==m1m1);
                }
                else if (cnns.size()==1) {
                    CHECK(cnns[0]==cnns2[0]);
                }
                else {
                    CHECK(cnns==cnns2);
                }
            }
        }
    }

}

TEST_CASE("emplace_before unit test"){
    //  0  : [ 4,  2,  3,  1,  5]
    //  1  : [ 7,  6,  2,  5,  0]
    //  2  : [ 7,  8,  3,  1,  0]
    //  3  : [ 9,  8,  4,  0,  2]
    //  4  : [ 9, 10,  3,  5,  0]
    //  5  : [ 6, 10,  4,  1,  0]
    //  6  : [11,  7, 10,  1,  5]
    //  7  : [11,  8,  6,  2,  1]
    //  8  : [11,  9,  7,  3,  2]
    //  9  : [11,  8, 10,  4,  3]
    //  10 : [11,  9,  6,  4,  5]
    //  11 : [ 9,  8,  7,  6, 10]

    SECTION("emplace_before at the beginning"){
        Triangulation<SPHERICAL_TRIANGULATION> trg(ICOSA_DATA, 0);
        trg.emplace_before(0, 4, 11);
        CHECK(trg.nodes().nn_ids(0)==std::vector<Index>{11, 4, 2, 3, 1, 5});
    }

    SECTION("emplace_before in the middle"){
        Triangulation<SPHERICAL_TRIANGULATION> trg(ICOSA_DATA, 0);
        trg.emplace_before(0, 2, 11);
        CHECK(trg.nodes().nn_ids(0)==std::vector<Index>{4, 11, 2, 3, 1, 5});
    }

    SECTION("emplace_before before the end"){
        Triangulation<SPHERICAL_TRIANGULATION> trg(ICOSA_DATA, 0);
        trg.emplace_before(0, 5, 11);
        CHECK(trg.nodes().nn_ids(0)==std::vector<Index>{4, 2, 3, 1, 11, 5});
    }

    SECTION("emplace_before nn not found"){
        Triangulation<SPHERICAL_TRIANGULATION> trg(ICOSA_DATA, 0);
        trg.emplace_before(0, 10, 11);
        CHECK(trg.nodes().nn_ids(0)==std::vector<Index>{4, 2, 3, 1, 5});
    }

}


TEST_CASE("unittest private member functions")
{
    SECTION("unittest all_nn_distance_vectors") {
        Triangulation<SPHERICAL_TRIANGULATION> star_triangulation(STAR_DATA, 0);
        star_triangulation.update_nn_distance_vectors(8);
        std::vector<vec3<Real>> all_nn_d = star_triangulation.nodes().nn_distances(8);//star_triangulation.all_nn_distance_vectors(8);
        auto v3 = vec3<Real>{1., 0., 0.};
        auto v1 = vec3<Real>{2., 0., 0.};
        auto v0 = vec3<Real>{0.5,0.,-1.};
        auto v6 = vec3<Real>{0.5, 0.866025_r, 0.0};
        auto v8 = vec3<Real>{1.5, 0.866025_r, 0.0};
        auto expected01 = v3 - v8;
        auto expected02 = v1 - v8;
        auto expected03 = v0 - v8;
        auto expected04 = v6 - v8;
        CHECK(all_nn_d[0]==expected01);
        CHECK(all_nn_d[1]==expected02);
        CHECK(all_nn_d[2]==expected03);
        CHECK(all_nn_d[3]==expected04);
    }

    SECTION("order_nn_ids: simple test on icosa data") {
        Triangulation<SPHERICAL_TRIANGULATION> icosa_triangulation_load(ICOSA_DATA, 0);
        CHECK(icosa_triangulation_load.order_nn_ids(0)==std::vector<Index>{3, 4, 5, 1, 2});
    }

    SECTION("orient_surface_of_a_sphere: simple test on icosa data") {
        Triangulation<SPHERICAL_TRIANGULATION> icosa_triangulation_load(ICOSA_DATA, 0);
        icosa_triangulation_load.orient_surface_of_a_sphere();

        CHECK(icosa_triangulation_load[0].nn_ids==std::vector<Index>{3, 4, 5, 1, 2});
        CHECK(icosa_triangulation_load[1].nn_ids==std::vector<Index>{6, 7, 2, 0, 5});
        CHECK(icosa_triangulation_load[3].nn_ids==std::vector<Index>{8, 9, 4, 0, 2});
        CHECK(icosa_triangulation_load[10].nn_ids==std::vector<Index>{9, 11, 6, 5, 4});

    }
}

TEST_CASE("unittest private static functions")
{
    SECTION("unittest cot_between_vectors for a few vectors") {
        auto li1 = vec3<Real>{1, 0, 0};
        auto li2 = vec3<Real>{0, -1, 0};
        Real cot_angle = Triangulation<SPHERICAL_TRIANGULATION>::cot_between_vectors(li1, li2);
        auto cot_angle_target = Catch::Approx(0).margin(0.001);
        CHECK(cot_angle==cot_angle_target);

        li1 = vec3<Real>{std::sqrt(3_r), 0, 0};
        li2 = vec3<Real>{std::sqrt(3_r), 0, 1};
        double cot_angle30 = Triangulation<SPHERICAL_TRIANGULATION>::cot_between_vectors(li1, li2);
        cot_angle_target = Catch::Approx(sqrt(3)).margin(0.001);
        CHECK(cot_angle30==cot_angle_target);

        li1 = vec3<Real>{1, std::sqrt(3_r), 0};
        li2 = vec3<Real>{1, 0, 0};
        double cot_angle60 = Triangulation<SPHERICAL_TRIANGULATION>::cot_between_vectors(li1, li2);
        cot_angle_target = Catch::Approx(1./sqrt(3)).margin(0.001);
        CHECK(cot_angle60==cot_angle_target);

        li1 = vec3<Real>{0, 1, 0};
        li2 = vec3<Real>{0, 1, 1};
        Real cot_angle45 = Triangulation<SPHERICAL_TRIANGULATION>::cot_between_vectors(li1, li2);
        cot_angle_target = Catch::Approx(1).margin(0.001);
        CHECK(cot_angle45==cot_angle_target);

        li1 = vec3<Real>{1, 0, 0};
        li2 = vec3<Real>{-1, 1, 0};
        Real cot_angleMin45 = Triangulation<SPHERICAL_TRIANGULATION>::cot_between_vectors(li1, li2);
        cot_angle_target = Catch::Approx(-1).margin(0.001);
        CHECK(cot_angleMin45==cot_angle_target);
    }
}