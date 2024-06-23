#ifndef FLIPPY_MONTECARLOUPDATER_HPP
#define FLIPPY_MONTECARLOUPDATER_HPP
/**
 * @file
 * @brief This file contains the MonteCarloUpdater class template. Together with Triangulation.hpp, this file contains flippy's most important high-level interfaces.
 */


#include "custom_concepts.hpp"
#include <random>
#include <variant>
#include "Nodes.hpp"
#include "Triangulation.hpp"

namespace fp {

/**
 * @brief A helper class for updating the triangulation, using
 * [Metropolis–Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm).
 * This is a high-level interface intended to implement a basic Monte-Carlo updating scheme, which should be sufficient
 * for a lot of simple situations.
 *
 * @tparam Real @RealStub
 * @tparam Index @IndexStub
 * @tparam EnergyFunctionParameters A user-defined struct type. An instance of this struct will be part of the input of the energy function.
 * Energy function provided during instantiation must have a third argument of type EnergyFunctionParameters const&.
 * @tparam RandomNumberEngine type of the random number engine that will be provided during the instantiation of the class.
 * For example [std::mt19937_64](https://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine).
 * @tparam triangulation_type One of the types specified by the TriangulationType enum.
 * This must match the type of triangulation provided during the class instantiation.
 */
template<typename EnergyFunctionParameters, typename RandomNumberEngine, TriangulationType triangulation_type>
class MonteCarloUpdater
{
private:
    static constexpr Real max_float = 3.40282347e+38f;
    fp::Triangulation<triangulation_type>& triangulation;
    EnergyFunctionParameters const& prms;
    typedef std::function<Real(fp::Node const&, fp::Triangulation<triangulation_type> const&, EnergyFunctionParameters const&, std::vector<Index> const&)> EnergyFunctionType;
    EnergyFunctionType energy_function;
    RandomNumberEngine& rng;
    std::uniform_real_distribution<Real> unif_distr_on_01;
    Real kBT_{1};
    Real min_bond_length_square{0.}, max_bond_length_square{max_float};
    unsigned long move_attempt{0}, bond_length_move_rejection{0},move_back{0};
    unsigned long flip_attempt{0}, bond_length_flip_rejection{0}, flip_back{0};

public:

    /**
     * @param triangulation_inp Reference to the triangulation that will be updated.
     * @param prms_inp The instance of the struct that contains the parameters of the system energy.
     * @param energy_function_inp A c++ function that represents the system energy. It can be evaluated for a given node for a Monte Carlo step.
     * @param rng_inp Random number engine
     * @param min_bond_length Minimal bond length between two triangulation nodes. More generally, a minimal distance two nodes of the triangulation are allowed to have.
     * If set to zero, the stability of the updater will suffer, and nonphysical shapes with self-intersecting triangulation will be common.
     * @param max_bond_length Maximal allowed length of a bond between two nodes of a triangulation.
     * If set too high, the stability of the updater will suffer, and nonphysical shapes with self-intersecting triangulation will be common.
     * Conversely, setting this variable too low will significantly reduce the chance of a successful bond flip.
     */
    MonteCarloUpdater(fp::Triangulation<triangulation_type>& triangulation_inp,
                      EnergyFunctionParameters const& prms_inp,
                      EnergyFunctionType energy_function_inp,
                      RandomNumberEngine& rng_inp, Real min_bond_length, Real max_bond_length)
    :triangulation(triangulation_inp), prms(prms_inp), energy_function(energy_function_inp), rng(rng_inp),
    unif_distr_on_01(std::uniform_real_distribution<Real>(0, 1)),
    min_bond_length_square(min_bond_length*min_bond_length), max_bond_length_square(max_bond_length*max_bond_length)
    {

    }

    //! Implementation of the Metropolis algorithm.
    /**
     * This function implements the [Metropolis algorithm](https://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm)
     * to evaluate whether the performed move should be rejected and undone.
     * @return `true` if the energy difference between the old and the new states is negative (i.e., the move costs energy)
     * and the random number is smaller than the Boltzmann Probability of accepting the move. `false` otherwise.
     */
    bool move_needs_undoing(Real e_diff)
    {
        if(kBT_>0){ //temperature can safely be put to 0, this will make the algorithm greedy
            return (e_diff<0) && (unif_distr_on_01(rng)>std::exp(e_diff/kBT_));
        }else{
            return (e_diff<0);
        }
    }

    //! Pre-update check to test that the update step will not result in an unphysical configuration.
    /**
     * Check if both the node does not overlap with any of its Verlet list neighbors or next neighbors and that none of
     * its connections to the next neighbors exceeds the maximal distance.
     * @param node @mcuNodeStub
     * @param displacement @mcuDisplacementStub
     * @return `true` if
     * new_next_neighbour_distances_are_between_min_and_max_length(fp::Node const&, fp::vec3<Real> const&)
     * and
     * new_verlet_neighbour_distances_are_between_min_and_max_length(fp::Node const&, fp::vec3<Real> const&)
     * conditions are both satisfied, `false` otherwise.
     */
    bool new_neighbour_distances_are_between_min_and_max_length(fp::Node const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        return (new_next_neighbour_distances_are_between_min_and_max_length(node, displacement)&&
            new_verlet_neighbour_distances_are_between_min_and_max_length(node, displacement));

    }

    //! Pre-update check to test that the update step will not result in an unphysical configuration.
    /**
     * Iterate through the next neighbor distances of a node which is a collection of
     * 3d vectors, i.e. nn_dist is of type fp::vec3. Then check if the displacement would make any of the
     * nn_distance vectors longer than the allowed max length.
     *
     * @param node @mcuNodeStub
     * @param displacement @mcuDisplacementStub
     * @return `true` if all next neighbor distances are between minimal and maximal allowed values,
     * provided during the instantiation of the MonteCarloUpdater class, `false` otherwise.
     */
    bool new_next_neighbour_distances_are_between_min_and_max_length(fp::Node const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        Real distance_square_new, distance_square_old;
        for (auto const& nn_dist: node.nn_distances){
            distance_square_new=(nn_dist - displacement).norm_square();
            distance_square_old=nn_dist.norm_square();
            if ((distance_square_new>max_bond_length_square) && (distance_square_old < max_bond_length_square)) {
                return false;
            }
            if ((distance_square_old>min_bond_length_square) && (distance_square_new<min_bond_length_square)) {
                return false;
            }
        }

        return true;
    }


    //! Pre-update check to test that the update step will not result in an unphysical configuration.
    /**
     * Iterate through the Verlet list neighbor distances of a node which is a collection of
     * 3d vectors. Then check if the displacement would make any of the Verlet list neighbors move closer than the minimal allowed distance
     * between two nodes. I.e., we check if the node would overlap with any of its Verlet list neighbors after the move.
     * The overlap distance or the minimal allowed distance is provided during the instantiation of the class.
     * @param node @mcuNodeStub
     * @param displacement @mcuDisplacementStub
     * @return `true` if the nodes did not overlap before but overlap now. `false` otherwise.
     */
    bool new_verlet_neighbour_distances_are_between_min_and_max_length(fp::Node const& node,
                                                                     fp::vec3<Real> const& displacement)

    {

        Real distance_square_new, distance_square_old;
        for (auto const& verlet_neighbour_id: node.verlet_list)
        {
            distance_square_new=(triangulation[verlet_neighbour_id].pos - node.pos - displacement).norm_square();
            distance_square_old=(triangulation[verlet_neighbour_id].pos - node.pos).norm_square();
            if ((distance_square_new<min_bond_length_square)&&(distance_square_old>min_bond_length_square)) { return false; }
        }
        return true;
    }

    //! Attempt a move Monte Carlo Step.
    /**
     * A move step is attempted on the specified node.
     * [Metropolis algorithm](https://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm) is used to evaluate whether the move should be accepted.
     * @param node @mcuNodeStub
     * @param displacement @mcuDisplacementStub
     * @return This function returns the calculated energy difference. If not energy difference was calculated due to bond_length constraint rejection, then an empty optional is returned.
     *
     */
    std::optional<Real> move_MC_updater(fp::Node const& node, fp::vec3<Real> const& displacement)
    {
        ++move_attempt;
        if (new_neighbour_distances_are_between_min_and_max_length(node, displacement)) {
            Real e_old = energy_function(node, triangulation, prms, node.nn_ids);
            triangulation.move_node(node.id, displacement);
            Real e_new = energy_function(node, triangulation, prms, node.nn_ids);
            Real e_diff = e_old - e_new;
            if (move_needs_undoing(e_diff)) {
                triangulation.move_node(node.id, -displacement);
                ++move_back;
            }
            return e_diff;
        }else{
            ++bond_length_move_rejection;
            return {};
        }
    }

    //! Attempt a flip Monte Carlo Step.
    /**
     * A flip step is attempted between a specified node and one of its randomly chosen next neighbors.
     * [Metropolis algorithm](https://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm) is used to evaluate whether the flip should be accepted.
     * @param node @mcuNodeStub
     * otherwise, the flippy will fail silently.
     * @note This function randomly chooses the next neighbor of the `node` and flips the edge between them.
     * If more precise control is required, i.e., if it is necessary to control exactly which edge needs to be flipped,
     * then the flip_MC_updater(fp::Node const& node, Index id_in_nn_ids) method can be used.
     */
    void flip_MC_updater(fp::Node const& node)
    {
      Index number_nn_ids = static_cast<Index>(node.nn_ids.size());
      Index nn_id = node.nn_ids[std::uniform_int_distribution<Index>(0, number_nn_ids-1)(rng)];
      flip_MC_updater(node, nn_id);
    }

    //! Attempt a flip Monte Carlo Step.
    /**
     * A flip step is attempted between a specified node and its specified next neighbor.
     * [Metropolis algorithm](https://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm) is used to evaluate whether the flip should be accepted.
     * @param node @mcuNodeStub
     * @param id_in_nn_ids This id must be part of the Node::nn_ids vector of the node provided in the first argument ,
     * otherwise, the flippy will fail silently.
     * @warning For performance reasons, this function does not check if `id_in_nn_ids` is really a `nn_id` of the `node` provided in the first argument.
     * The user is required to guarantee this fact. Otherwise, flippy will fail unpredictably.
     * flip_MC_updater(fp::Node const&) is the safer method if the user does not care exactly which bond is flipped.
     */
    void flip_MC_updater(fp::Node const& node, Index id_in_nn_ids)
    {
        ++flip_attempt;
        Neighbors cns = triangulation.previous_and_next_neighbour_global_ids(node.id, id_in_nn_ids);
        //This one line leads to an allocation each MC step for each node, which can easily number in billions
        // for even small simulations. This can not be deploied to production! For even single flip per move, this
        // leads to a 10% regression in banchmarks.
        const auto changed_neighborhood = std::vector{id_in_nn_ids, cns.j_m_1, cns.j_p_1};
        Real e_old = energy_function(node, triangulation, prms, changed_neighborhood);
        BondFlipData bfd = triangulation.flip_bond(node.id, id_in_nn_ids, min_bond_length_square, max_bond_length_square);
        if (bfd.flipped) {
            Real e_new = energy_function(node, triangulation, prms,  changed_neighborhood);
            Real e_diff = e_old - e_new;
            if (move_needs_undoing(e_diff)) {
                triangulation.unflip_bond(node.id, id_in_nn_ids, bfd);
                ++flip_back;

            }
        }else{++bond_length_flip_rejection;}
    }

    //! Reset the temperature of the Monte Carlo updater, at which the Boltzmann weights are evaluated.
    void reset_kBT(Real kBT){
        /**
         * @param kBT input value of new temperature.
         * The internal private data member `kBT_` that contains the current value of temperature will be overwritten
         * with this input.
         */
        kBT_=kBT;
    }

    //! @getterFunctionStub
    [[nodiscard]] Real kBT(){
    /**
     * Monte Carlo Updater requires a kBT value in the calculation of the rejection probability of a random move.
     * @return the state of the kBT value of the updater.
     * @see move_needs_undoing()
     */
        return kBT_;
    }
    //! @getterFunctionStub
    [[nodiscard]] unsigned long move_attempt_count() const {
    /**
     * Every time a move is attempted, a private internal state variable `move_attempt` is incremented by move_MC_updater().
     * This variable can be used for diagnostics or statistical tracking, but its state does not impact the function of the updater.
     * @return current state of `move_attempt`.
     * @see move_needs_undoing()
     */
        return move_attempt;
    }
    //! @getterFunctionStub
    [[nodiscard]] unsigned long bond_length_move_rejection_count() const {
    /**
     * Moves that cause nodes to overlap with their Verlet list neighbors or move them too far away from any of their next neighbors are rejected. Every time such rejection happens, a private internal state variable `bond_length_move_rejection` is incremented by move_MC_updater().
     * The specifics of this rejection criteria are calculated in the function new_neighbour_distances_are_between_min_and_max_length(fp::Node const& node, fp::vec3<Real> const& displacement)
     * Every time a node move would lead that node to have a bond with one of its next neighbors, which is longer than a specified maximal length (max_bond_length()), the move will be rejected.
     * This variable can be used for diagnostics or statistical tracking, but its state does not impact the function of the updater.
     * @return current state of `bond_length_move_rejection`.
     * @see move_needs_undoing()
     */
        return bond_length_move_rejection;
    }
    //! @getterFunctionStub
    [[nodiscard]] unsigned long move_back_count() const {
    /**
     * Every time a move is rejected because the energy requirement was not satisfied, a private internal state variable `move_back` is incremented by move_MC_updater().
     * This variable does not track the rejections resulting from bond length restriction violations.
     * This variable can be used for diagnostics or statistical tracking, but its state does not impact the function of the updater.
     * @return current state of `move_back`.
     * @see move_needs_undoing() bond_length_move_rejection_count()
     */
        return move_back;
    }
    //!@getterFunctionStub
    [[nodiscard]] unsigned long flip_attempt_count() const {
    /**
     * Every time a flip is attempted, a private internal state variable `flip_attempt` is incremented by flip_MC_updater() and flip_MC_updater(fp::Node const& node, Index index_in_nn_ids).
     * This variable can be used for diagnostics or statistical tracking, but its state does not impact the function of the updater.
     * @return current state of `flip_attempt`.
     */
        return flip_attempt;
    }
    //! @getterFunctionStub
    [[nodiscard]] unsigned long bond_length_flip_rejection_count() const {
    /**
     * If a flip would turn a valid bond into a bond that is too long, the flip is rejected, a private internal state variable `bond_length_flip_rejection` is incremented by flip_MC_updater() and flip_MC_updater(fp::Node const& node, Index index_in_nn_ids).
     * The rejection of flips is handled by the Triangulation class itself and reported through the BondFlipData to the flip functions of the updater.
     * This variable can be used for diagnostics or statistical tracking, but its state does not impact the function of the updater.
     * @return current state of `bond_length_flip_rejection`.
     * @see Triangulation::flip_bond(Index, Index, Real, Real) Triangulation::unflip_bond(Index, Index, BondFlipData const&)
     */
        return bond_length_flip_rejection;
    }
    //! @getterFunctionStub
    [[nodiscard]] unsigned long flip_back_count() const {
    /**
     * Every time a flip is rejected because the energy requirement was not satisfied, a private internal state variable `flip_back` is incremented by flip_MC_updater() and flip_MC_updater(fp::Node const& node, Index index_in_nn_ids).
     * The rejection of flips is handled by the Triangulation class itself and reported through the BondFlipData to the flip functions of the updater.
     * This variable does not track the rejections resulting from bond length restriction violations.
     * This variable can be used for diagnostics or statistical tracking, but its state does not impact the function of the updater.
     * @return current state of `flip_back`.
     * @see Triangulation:: unflip_bond(Index node_id, Index nn_id, BondFlipData const& common_nns) bond_length_flip_rejection_count()
     */
        return flip_back;
    }

};
}
#endif //FLIPPY_MONTECARLOUPDATER_HPP
