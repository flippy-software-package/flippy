#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "flippy.hpp"
#include "test_helper.h"

TEMPLATE_TEST_CASE("util vol and area tests", "", float, double ){

    using REAL = TestType;
    SECTION("spherical volume"){
        REAL R = 17.1f;
        CHECK_THAT(fp::sphere_vol<REAL>(R), 20944.834858665166_ish);
    }

    SECTION("spherical area"){
        REAL R = 11.1f;
            CHECK_THAT(fp::sphere_area<REAL>(R), 1548.3025233951935_ish);
    }

}

TEMPLATE_TEST_CASE("linear adaptation function test", "", float, double){
    using REAL = TestType;
    REAL t_init = 1.f;
    REAL t_fin = 3.f;
    REAL v_init = 2.f;
    REAL v_fin = 8.f;
    SECTION("linear adaptation val_init plateau"){
        CHECK_THAT(fp::linear_adaptation<REAL>(0.f, t_init, t_fin, v_init, v_fin), ish(v_init));
    }

    SECTION("linear adaptation val_fin plateau"){
        CHECK_THAT(fp::linear_adaptation<REAL>(2.f*t_fin, t_init, t_fin, v_init, v_fin), ish(v_fin));
    }

    SECTION("linear adaptation middle value"){
        REAL t_middle = t_init + 0.5f*(t_fin - t_init);
        REAL v_middle = v_init + 0.5f*(v_fin - v_init);
        CHECK_THAT(fp::linear_adaptation<REAL>(t_middle, t_init, t_fin, v_init, v_fin), ish(v_middle));
    }

    SECTION("decrease also works"){
        REAL v_fin2 = 1.f;
        REAL v_init2 = 5.f;
        REAL t_middle = t_init + 0.5f*(t_fin - t_init);
        REAL v_middle = v_init2 + 0.5f*(v_fin2 - v_init2);
        CHECK_THAT(fp::linear_adaptation<REAL>(t_middle, t_init, t_fin, v_init2, v_fin2), ish(v_middle));
    }


    SECTION("calculation of needed steps for fixed speed adaptation 0"){
        REAL t0 = 0.f;
        REAL val_init = 1.f;
        REAL val_fin = 0.5f;
        REAL delta_val = -0.1f;
        auto t_steps = fp::steps_for_fixed_speed_adaptation<REAL, unsigned int>(t0, val_init, val_fin, delta_val);
        unsigned int expected = 5;

        CHECK(t_steps == expected);
    }

    SECTION("calculation of needed steps for fixed speed adaptation 1"){
        REAL t0 = 1.f;
        REAL val_init = 1.f;
        REAL val_fin = 0.5f;
        REAL delta_val = -0.1f;
        auto t_steps = fp::steps_for_fixed_speed_adaptation<REAL, unsigned int>(t0, val_init, val_fin, delta_val);
        unsigned int expected = 6;

        CHECK(t_steps == expected);
    }

    SECTION("calculation of needed steps for fixed speed adaptation, calculate one more step if delta_cal is not a clean divisor of val_interval"){
        REAL t0 = 0.f;
        REAL val_init = 1.f;
        REAL val_fin = 0.45f;
        REAL delta_val = -0.1f;
        auto t_steps = fp::steps_for_fixed_speed_adaptation<REAL, unsigned int>(t0, val_init, val_fin, delta_val);
        unsigned int expected = 6;

        CHECK(t_steps == expected);
    }

}
