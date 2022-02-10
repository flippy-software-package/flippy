#ifndef FLIPPY_MONTECARLOUPDATER_HPP
#define FLIPPY_MONTECARLOUPDATER_HPP
#include <concepts>
#include <random>
#include "Nodes.hpp"
#include "Triangulation.hpp"

namespace fp {

template<std::floating_point Real, std::integral Index, typename EnergyFunctionParameters, typename RandomNumberEngine, TriangulationType triangulation_type>
class MonteCarloUpdater
{
private:
    Real e_old{}, e_new{}, e_diff{};
    fp::Triangulation<Real, Index, triangulation_type>& triangulation;
    EnergyFunctionParameters const& prms;
    std::function<Real(fp::Node<Real, Index> const&, fp::Triangulation<Real, Index, triangulation_type> const&, EnergyFunctionParameters const&)> energy_function;
    RandomNumberEngine& rng;
    std::uniform_real_distribution<Real> unif_distr_on_01;
    Real min_bond_length_square, max_bond_length_square;
    Index move_attempt{0}, bond_length_move_rejection{0},move_back{0};
    Index flip_attempt{0}, bond_length_flip_rejection{0}, flip_back{0};

public:
    MonteCarloUpdater()=default;
    MonteCarloUpdater(fp::Triangulation<Real, Index, triangulation_type>& triangulation_inp,
                      EnergyFunctionParameters const& prms_inp,
                      std::function<Real(fp::Node<Real, Index> const&, fp::Triangulation<Real, Index, triangulation_type> const&, EnergyFunctionParameters const&)> energy_function_inp,
                      RandomNumberEngine& rng_inp, Real min_bond_length, Real max_bond_length)
    :triangulation(triangulation_inp), prms(prms_inp), energy_function(energy_function_inp), rng(rng_inp),
    unif_distr_on_01(std::uniform_real_distribution<Real>(0, 1)),
    min_bond_length_square(min_bond_length*min_bond_length), max_bond_length_square(max_bond_length*max_bond_length)
    {}

    bool move_needs_undoing()
    {
        e_diff = e_old - e_new;
        return (e_diff<0) && (unif_distr_on_01(rng)>std::exp(e_diff));
    }

    bool new_neighbour_distances_are_between_min_and_max_length(fp::Node<Real, Index> const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        return (new_next_neighbour_distances_are_between_min_and_max_length(node, displacement)&&
            new_verlet_neighbour_distances_are_between_min_and_max_length(node, displacement));

    }

    bool new_next_neighbour_distances_are_between_min_and_max_length(fp::Node<Real, Index> const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        /**
         * iterate through the next neighbour distances of a node which is a collection of
         * 3d vectors, i.e. nn_dist is of type fp::vec3. Then check if the displacement would make any of the
         * nn_distance vectors longer than the allowed max length.
         */
        Real distance_square;
        for (auto const& nn_dist: node.nn_distances)
        {
            distance_square=(nn_dist - displacement).norm_square();
            if ((distance_square>max_bond_length_square)||(distance_square<min_bond_length_square)) { return false; }
        }
        return true;
    }

    bool new_verlet_neighbour_distances_are_between_min_and_max_length(fp::Node<Real, Index> const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        /**
         * iterate through the next neighbour distances of a node which is a collection of
         * 3d vectors, i.e. nn_dist is of type fp::vec3. Then check if the displacement would make any of the
         * nn_distance vectors longer than the allowed max length.
         */
        Real distance_square;
        for (auto const& verlet_neighbour_id: node.verlet_list)
        {
            distance_square=(triangulation[verlet_neighbour_id].pos - displacement).norm_square();
            if (distance_square<min_bond_length_square) { return false; }
        }
        return true;
    }

    void move_MC_updater(fp::Node<Real, Index> const& node, fp::vec3<Real> const& displacement)
    {
        ++move_attempt;
        if (new_neighbour_distances_are_between_min_and_max_length(node, displacement)) {
            e_old = energy_function(node, triangulation, prms);
            triangulation.move_node(node.id, displacement);
            e_new = energy_function(node, triangulation, prms);
            if (move_needs_undoing()) {triangulation.move_node(node.id, -displacement); ++move_back;}
        }else{++bond_length_move_rejection;}
    }

    void flip_MC_updater(fp::Node<Real, Index> const& node)
    {
        ++flip_attempt;
        e_old = energy_function(node, triangulation, prms);
        Index number_nn_ids = node.nn_ids.size();
        Index nn_id = node.nn_ids[std::uniform_int_distribution<Index>(0, number_nn_ids-1)(rng)];
        auto bfd = triangulation.flip_bond(node.id, nn_id, min_bond_length_square, max_bond_length_square);
        if (bfd.flipped) {
            e_new = energy_function(node, triangulation, prms);
            if (move_needs_undoing()) { triangulation.unflip_bond(node.id, nn_id, bfd); ++flip_back;}
        }else{++bond_length_flip_rejection;}
    }

};
}
#endif //FLIPPY_MONTECARLOUPDATER_HPP
