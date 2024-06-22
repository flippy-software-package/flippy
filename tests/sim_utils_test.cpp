#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "flippy.hpp"
#include "test_helper.h"

TEST_CASE("util vol and area tests"){

    SECTION("spherical volume"){
        fp::Real R = 17.1f;
        CHECK_THAT(fp::sphere_vol(R), 20944.834858665166_ish);
    }

    SECTION("spherical area"){
        fp::Real R = 11.1f;
        CHECK_THAT(fp::sphere_area(R), 1548.3025233951935_ish);
    }

}

TEST_CASE("linear adaptation function test"){
    fp::Real t_init = 1.f;
    fp::Real t_fin = 3.f;
    fp::Real v_init = 2.f;
    fp::Real v_fin = 8.f;
    SECTION("linear adaptation val_init plateau"){
        CHECK_THAT(fp::linear_adaptation(0.f, t_init, t_fin, v_init, v_fin), ish(v_init));
    }

    SECTION("linear adaptation val_fin plateau"){
        CHECK_THAT(fp::linear_adaptation(2.f*t_fin, t_init, t_fin, v_init, v_fin), ish(v_fin));
    }

    SECTION("linear adaptation middle value"){
        fp::Real t_middle = t_init + 0.5f*(t_fin - t_init);
        fp::Real v_middle = v_init + 0.5f*(v_fin - v_init);
        CHECK_THAT(fp::linear_adaptation(t_middle, t_init, t_fin, v_init, v_fin), ish(v_middle));
    }

    SECTION("decrease also works"){
        fp::Real v_fin2 = 1.f;
        fp::Real v_init2 = 5.f;
        fp::Real t_middle = t_init + 0.5f*(t_fin - t_init);
        fp::Real v_middle = v_init2 + 0.5f*(v_fin2 - v_init2);
        CHECK_THAT(fp::linear_adaptation(t_middle, t_init, t_fin, v_init2, v_fin2), ish(v_middle));
    }


    SECTION("calculation of needed steps for fixed speed adaptation 0"){
        fp::Real t0 = 0.f;
        fp::Real val_init = 1.f;
        fp::Real val_fin = 0.5f;
        fp::Real delta_val = -0.1f;
        auto t_steps = fp::steps_for_fixed_speed_adaptation(t0, val_init, val_fin, delta_val);
        unsigned int expected = 5;

        CHECK(t_steps == expected);
    }

    SECTION("calculation of needed steps for fixed speed adaptation 1"){
        fp::Real t0 = 1.f;
        fp::Real val_init = 1.f;
        fp::Real val_fin = 0.5f;
        fp::Real delta_val = -0.1f;
        auto t_steps = fp::steps_for_fixed_speed_adaptation(t0, val_init, val_fin, delta_val);
        unsigned int expected = 6;

        CHECK(t_steps == expected);
    }

    SECTION("calculation of needed steps for fixed speed adaptation, calculate one more step if delta_cal is not a clean divisor of val_interval"){
        fp::Real t0 = 0.f;
        fp::Real val_init = 1.f;
        fp::Real val_fin = 0.45f;
        fp::Real delta_val = -0.1f;
        auto t_steps = fp::steps_for_fixed_speed_adaptation(t0, val_init, val_fin, delta_val);
        unsigned int expected = 6;

        CHECK(t_steps == expected);
    }

}
