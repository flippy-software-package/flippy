#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <iostream>
#define TESTING_FLIPPY_TRIANGULATION_ndh6jclc0qnp274b = 1
#include "flippy.hpp"
using namespace fp;
const double EPSILON = 1e-9;

static fp::Json ICOSA_DATA(){
        return R"({
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
}

void rescale_triangulation(Real R, Triangulation& tr)
{
    tr.R_initial=R;
    tr.scale_all_nodes_to_R_init();
    tr.orient_surface_of_a_sphere();
    tr.initiate_distance_vectors();
    tr.make_global_geometry();
    tr.make_verlet_list();
}
TEST_CASE("Icosa geometry check test")
{
    Real R = 20._r;
    Triangulation icosahedron(ICOSA_DATA(),  0);
    rescale_triangulation(R, icosahedron);

    Real sin_ = std::sin(2_r*PI/5._r);
    Real sin_2 = sin_*sin_;
    Real A_SUB = R*R/(4._r*std::tan(PI/3_r)*sin_2 );
    auto A_NODE = 5*A_SUB;
    auto A_NODE_target = Catch::Approx(A_NODE).margin(EPSILON);
    Real V_SUB = (R/3_r)*A_SUB*std::sqrt(sin_2 - 1_r/3._r)/sin_;
    auto V_NODE = 5*V_SUB;
    auto V_NODE_target = Catch::Approx(V_NODE).margin(EPSILON);
    auto K_SQUARE_NODE = 0.5*std::pow(2_r/R, 2_r)*A_NODE;
    auto K_SQUARE_NODE_target = Catch::Approx(K_SQUARE_NODE).margin(EPSILON);
    for (auto const& node : icosahedron.nodes().data) {
        CHECK(A_NODE_target==node.area);
        CHECK(node.volume==V_NODE_target);
        CHECK(node_unit_bending_energy(node)==K_SQUARE_NODE_target);
    }
}

TEST_CASE("Test Geometry Container")
{

    SECTION("uninitiated") {
        struct GeometryWrapper {
          Geometry geo;
        };
        GeometryWrapper gw;
        CHECK(gw.geo.area==0.);
        CHECK(gw.geo.volume==0.);
    }

    SECTION("default initiator") {
        Geometry geo{};
        CHECK(geo.area==0.);
        CHECK(geo.volume==0.);
    }

}

TEST_CASE("Sphere geometry test")
{
    Real R = 1000.;
    auto all_data = json_read("../../tests/init_files/egg_init.json");
    Triangulation sphere(all_data["nodes"], 0);
    rescale_triangulation(R, sphere);
    auto ACCEPTABLE_ERROR = 0.01;

    Real A_SPHERE = 4_r*PI*R*R;
    auto A_SPHERE_target = Catch::Approx(A_SPHERE).epsilon(ACCEPTABLE_ERROR);

    Real V_SPHERE = A_SPHERE*R/3_r;
    auto V_SPHERE_target = Catch::Approx(V_SPHERE).epsilon(ACCEPTABLE_ERROR);

    auto UNIT_BENDING_ENERGY_SPHERE = 8_r*PI;
    auto UNIT_BENDING_ENERGY_SPHERE_target = Catch::Approx(UNIT_BENDING_ENERGY_SPHERE).epsilon(ACCEPTABLE_ERROR);

    Geometry lg{};
    Real aggregated_energy = 0;
    for (auto const& node : sphere.nodes()) {
        lg += Geometry(node);
        aggregated_energy += node_unit_bending_energy(node);
    }

    SECTION("externally calculated global geometry") {
        CHECK(lg.area==A_SPHERE_target);
        CHECK(lg.volume==V_SPHERE_target);
        CHECK(aggregated_energy==UNIT_BENDING_ENERGY_SPHERE_target);
    }

    SECTION("internally calculated global geometry") {
        CHECK(sphere.global_geometry().area==A_SPHERE_target);
        CHECK(sphere.global_geometry().volume==V_SPHERE_target);
    }

    vec3<Real> displacement{1200., 2329., -12.901_r};
    sphere.translate_all_nodes(displacement);

    Geometry lg_translated{};
    Real aggregated_energy_translated = 0;
    for (auto const& node : sphere.nodes()) {
        lg_translated += Geometry(node);
        aggregated_energy_translated += node_unit_bending_energy(node);
    }

    SECTION("checking translational invariance") {
        CHECK(lg_translated.area==A_SPHERE_target);
        CHECK(lg_translated.volume==V_SPHERE_target);
        CHECK(aggregated_energy_translated==UNIT_BENDING_ENERGY_SPHERE_target);
    }

    SECTION("checking translational invariance") {
        CHECK(sphere.global_geometry().area==A_SPHERE_target);
        CHECK(sphere.global_geometry().volume==V_SPHERE_target);
    }
}

TEST_CASE("Ellipse geometry test")
{
    Real R = 10.;
    auto all_data = json_read("../../tests/init_files/egg_init.json");
    double x_stretch = 1.4;
    Triangulation ellipse(all_data["nodes"], 0);
    rescale_triangulation(R, ellipse);
    ellipse.scale_node_coordinates(static_cast<Real>(x_stretch));

    auto ACCEPTABLE_ERROR = 0.01;

    Real V_Ellipse = static_cast<Real>(x_stretch*4.*M_PI*R*R*R/3.);
    auto V_Ellipse_target = Catch::Approx(V_Ellipse).epsilon(ACCEPTABLE_ERROR);
    double e = std::sqrt(1. - 1./(x_stretch*x_stretch));
    auto A_Ellipse = static_cast<Real>(2.*M_PI*R*R*(1 + (x_stretch/e)*std::asin(e)));
    auto A_Ellipse_target = Catch::Approx(A_Ellipse).epsilon(ACCEPTABLE_ERROR);

//    json_dump("../../../../data/ellipse_egg", ellipse.make_egg_data());

    SECTION("internally calculated global geometry") {
        CHECK(ellipse.global_geometry().volume==V_Ellipse_target);
        CHECK(ellipse.global_geometry().area==A_Ellipse_target);
    }

    vec3<Real> displacement{1200., 2329., -12.901_r};
    ellipse.translate_all_nodes(displacement);


    SECTION("checking translational invariance1") {
        CHECK(ellipse.global_geometry().volume==V_Ellipse_target);
        CHECK(ellipse.global_geometry().area==A_Ellipse_target);
    }

    vec3<Real> displacement2{0., 9., -1};
    ellipse.translate_all_nodes(displacement2 - displacement);


    SECTION("checking translational invariance2") {
        CHECK(ellipse.global_geometry().volume==V_Ellipse_target);
        CHECK(ellipse.global_geometry().area==A_Ellipse_target);
    }
}

TEST_CASE("Ellipse geometry test for triangulator mesh")
{
    Real R = 2.;
    Real x_stretch = 1.3_r;
    Triangulation ellipse(15, R, 0);
    ellipse.scale_node_coordinates(x_stretch);

    auto ACCEPTABLE_ERROR = 0.02;

    Real V_Ellipse = x_stretch*4_r*PI*R*R*R/3_r;
    auto V_Ellipse_target = Catch::Approx(V_Ellipse).epsilon(ACCEPTABLE_ERROR);
    auto e = std::sqrt(1_r - 1_r/(x_stretch*x_stretch));
    Real A_Ellipse = 2*PI*R*R*(1 + (x_stretch/e)*asin(e));
    auto A_Ellipse_target = Catch::Approx(A_Ellipse).epsilon(ACCEPTABLE_ERROR);

    SECTION("internally calculated global geometry") {
        CHECK(ellipse.global_geometry().volume==V_Ellipse_target);
        CHECK(ellipse.global_geometry().area==A_Ellipse_target);
    }

    vec3<Real> displacement{1200., 2329., -12.901_r};
    ellipse.translate_all_nodes(displacement);


    SECTION("checking translational invariance1") {
        CHECK(ellipse.global_geometry().area==A_Ellipse_target);
        CHECK(ellipse.global_geometry().volume==V_Ellipse_target);
    }

    vec3<Real> displacement2{0., 9., -1};
    ellipse.translate_all_nodes(displacement2 - displacement);


    SECTION("checking translational invariance") {
        CHECK(ellipse.global_geometry().volume==V_Ellipse_target);
        CHECK(ellipse.global_geometry().area==A_Ellipse_target);
    }
}

static fp::Json CUBE_DATA(){
        return R"({
	  "0":	{"nn_ids": [3,2,1,4],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,0]},
	  "1":  {"nn_ids": [0,2,6,5,4], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,0]},
	  "2":  {"nn_ids": [1,0,3,7,6], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [2,2,0]},
	  "3":  {"nn_ids": [0,4,7,2],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [2,0,0]},
	  "4":  {"nn_ids": [3,0,1,5,7], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,-2]},
	  "5":  {"nn_ids": [6,7,4,1],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,-2]},
	  "6":  {"nn_ids": [1,2,7,5],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [2,2,-2]},
	  "7":  {"nn_ids": [6,2,3,4,5], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [2,0,-2]}
  })"_json;
}

static fp::Json BRICK_DATA()
{
    return R"({
	  "0":	{"nn_ids": [3,2,1,4],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,0]},
	  "1":  {"nn_ids": [0,2,6,5,4], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,0]},
	  "2":  {"nn_ids": [1,0,3,7,6], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [3,2,0]},
	  "3":  {"nn_ids": [0,4,7,2],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [3,0,0]},
	  "4":  {"nn_ids": [3,0,1,5,7], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,-2]},
	  "5":  {"nn_ids": [6,7,4,1],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,-2]},
	  "6":  {"nn_ids": [1,2,7,5],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [3,2,-2]},
	  "7":  {"nn_ids": [6,2,3,4,5], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [3,0,-2]}
  })"_json;
}

TEST_CASE("cube geometry test")
{
    Triangulation cube(CUBE_DATA(), 0);
    Real ACCEPTABLE_ERROR = 0.01_r;
    auto A_square = Catch::Approx(24.).margin(ACCEPTABLE_ERROR);
    auto V_square = Catch::Approx(8).margin(ACCEPTABLE_ERROR);
    SECTION("cube geometry") {
        CHECK(cube.global_geometry().area==A_square);
        CHECK(cube.global_geometry().volume==V_square);
    }

    cube.translate_all_nodes({10, -187.1_r, -9.1_r});
    SECTION("cube geometry after shift") {
        CHECK(cube.global_geometry().area==A_square);
        CHECK(cube.global_geometry().volume==V_square);
    }

}

TEST_CASE("Brick geometry test")
{
    Triangulation brick(BRICK_DATA(), 0);
    Real ACCEPTABLE_ERROR = 0.01_r;
    auto A_square = Catch::Approx(32.).margin(ACCEPTABLE_ERROR);
    auto V_square = Catch::Approx(12).margin(ACCEPTABLE_ERROR);
    SECTION("square geometry") {
        CHECK(brick.global_geometry().area==A_square);
        CHECK(brick.global_geometry().volume==V_square);
    }

    brick.translate_all_nodes({10_r, -187.1_r, -9.1_r});
    SECTION("brick geometry after shift") {
        CHECK(brick.global_geometry().area==A_square);
        CHECK(brick.global_geometry().volume==V_square);
    }

}

static fp::Json HYPERBRICK_DATA()
{
    return R"({
	  "0":	{"nn_ids": [3,2,1,4],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,0]},
	  "1":  {"nn_ids": [0,2,6,5,4], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,0]},
	  "2":  {"nn_ids": [1,0,3,7,6], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [300,2,0]},
	  "3":  {"nn_ids": [0,4,7,2],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [300,0,0]},
	  "4":  {"nn_ids": [3,0,1,5,7], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,-2]},
	  "5":  {"nn_ids": [6,7,4,1],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,-2]},
	  "6":  {"nn_ids": [1,2,7,5],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [300,2,-2]},
	  "7":  {"nn_ids": [6,2,3,4,5], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [300,0,-2]}
  })"_json;
}

static fp::Json HYPER_SQUEEZED_BRICK_DATA(){
    return R"({
	  "0":	{"nn_ids": [3,2,1,4],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,0]},
	  "1":  {"nn_ids": [0,2,6,5,4], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,0]},
	  "2":  {"nn_ids": [1,0,3,7,6], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.03,2,0]},
	  "3":  {"nn_ids": [0,4,7,2],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.03,0,0]},
	  "4":  {"nn_ids": [3,0,1,5,7], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,0,-2]},
	  "5":  {"nn_ids": [6,7,4,1],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0,2,-2]},
	  "6":  {"nn_ids": [1,2,7,5],   "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.03,2,-2]},
	  "7":  {"nn_ids": [6,2,3,4,5], "verlet_list": [], "curvature_vec": [0,0,0], "area": 0, "volume": 0, "unit_bending_energy": 0, "pos": [0.03,0,-2]}
  })"_json;
}

TEST_CASE("Hyperstretch geometry test")
{
    Triangulation brick(HYPERBRICK_DATA(), 0);
    Real ACCEPTABLE_ERROR = 0.01_r;
    auto A_square = Catch::Approx(2408.).margin(ACCEPTABLE_ERROR);
    auto V_square = Catch::Approx(1200).margin(ACCEPTABLE_ERROR);
    SECTION("square geometry") {
        CHECK(brick.global_geometry().area==A_square);
        CHECK(brick.global_geometry().volume==V_square);
    }

    brick.translate_all_nodes({10_r, -187.1_r, -9.1_r});
    SECTION("brick geometry after shift") {
        CHECK(brick.global_geometry().area==A_square);
        CHECK(brick.global_geometry().volume==V_square);
    }
}

TEST_CASE("Hyperqueeze geometry test")
{
    Triangulation brick(HYPER_SQUEEZED_BRICK_DATA(), 0);
    Real ACCEPTABLE_ERROR = 0.01_r;
    auto A_square = Catch::Approx(8.24).margin(ACCEPTABLE_ERROR);
    auto V_square = Catch::Approx(0.12).margin(ACCEPTABLE_ERROR);
    SECTION("square geometry") {
        CHECK(brick.global_geometry().area==A_square);
        CHECK(brick.global_geometry().volume==V_square);
    }

    brick.translate_all_nodes({0_r, -17.1_r, 90.1_r});
    SECTION("brick geometry after shift") {
        CHECK(brick.global_geometry().area==A_square);
        CHECK(brick.global_geometry().volume==V_square);
    }
}