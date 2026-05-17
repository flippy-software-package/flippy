#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "catch2/generators/catch_generators_all.hpp"
#include "flippy.hpp"

using namespace fp;

bool euler_number_is_2(const Triangulation &triangulation) {
    /**
     * General euler formula for flat polyhedra: `V - E + F = 2`
     * Special formula for triangulations is derived through `F = (2/3)E`:
     *
     * `E = 3V - 6`
     *
     */
    Index V = triangulation.size();
    Index E = 0;

    for (const auto &node : triangulation.nodes()) {
        E += node.nn_ids.size();
    }

    E = E / 2;

    return E == 3 * V - 6;
}

SCENARIO("Random bond flips do not destroy the topology") {
    GIVEN("A a spherical triangulation") {
        constexpr Index subriang_level   = 10;
        constexpr auto  verlet_radius    = 0_r;
        constexpr auto  r_init           = 13_r;
        constexpr Index test_repetitions = 10;

        auto tr = Triangulation::make_spherical_triangulation(subriang_level, r_init, verlet_radius);
        auto rg = Catch::Generators::RandomIntegerGenerator<Index>(0, tr.size(), Catch::Generators::Detail::getSeed());

        THEN("Euler Number is 2") { CHECK(euler_number_is_2(tr)); }

        WHEN("Bonds are flipped randomly") {
            Index rn = rg.get();
            rg.get();

            for (Index _ = 0; _ < test_repetitions; ++_) {
                for (const auto &node : tr.nodes()) {
                    auto  n_neighbors = static_cast<Index>(node.nn_ids.size());
                    Index rand_id     = rn % n_neighbors;
                    tr.flip_bond(node.id, node.nn_ids[rand_id], 0, r_init * r_init);
                }
            }

            THEN("The the euler characteristic of the triangulation does not change") { CHECK(euler_number_is_2(tr)); }
        }
    }
}
