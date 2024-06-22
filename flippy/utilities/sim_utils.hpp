#ifndef FLIPPY_SIM_UTILS_H
#define FLIPPY_SIM_UTILS_H
#include <custom_concepts.hpp>
#include <random>
#include <optional>
#include <Nodes.hpp>
#include <vec3.hpp>

namespace fp{

    [[maybe_unused]] [[nodiscard]] static Real sphere_vol(Real R){return (4._r/3._r)* PI *R*R*R;}
    [[maybe_unused]] [[nodiscard]] static Real sphere_area(Real R){return 4._r * PI *R*R;}

    [[maybe_unused]] [[nodiscard]] static Real linear_adaptation(Real x, Real x_init, Real x_fin,
            Real val_init, Real val_fin)
     {
        if (x <= x_init) { return val_init; }
        if (x > x_fin) { return val_fin; }
//        if (std::abs(x_fin-x_init)<0.001){return val_init;}
        return val_init + (val_fin-val_init)/(x_fin-x_init)*(x-x_init);
    }

    [[maybe_unused]] [[nodiscard]] static Index steps_for_fixed_speed_adaptation(Real x_init, Real val_init, Real val_fin, Real delta_val)
    {
        return static_cast<Index>(std::ceil((val_fin-val_init)/delta_val + x_init));
    }

    //! N_node=12+30*n+20*n*(n-1)/2 where n is the same as sub_triangulation_iteration_count
    [[maybe_unused]] [[nodiscard]] static Real min_radius_with_non_overlapping_beads(Real min_allowed_distance_between_bead_centers, Index sub_triangulation_iteration_count){
        Real l_min = min_allowed_distance_between_bead_centers;
        Real n = static_cast<Real>(sub_triangulation_iteration_count);

        return l_min / (2_r * std::sin(std::asin(1_r / (2_r * std::sin(2_r * PI / 5_r ))) / (n + 1_r)));
    }

    class DynamicDisplacementUpdater{
        Real d0_;
        Real delta_d_;
        Real p_target_;
        Real p_accum_{static_cast<Real>(0.f)};
        Real p_count_{static_cast<Real>(1.f)};

        Real prob(std::optional<Real> e_diff, Real kBT){ return (e_diff.has_value()?std::min(std::exp(e_diff.value()/kBT),1_r):0_r); }

    public:
        DynamicDisplacementUpdater(Real initial_displacement, Real displacement_adaptation_step_size, Real p_target)
        : d0_(initial_displacement), delta_d_(displacement_adaptation_step_size), p_target_(p_target) { }

        Real new_displacement_magnitude(){
            Real p = p_accum_/p_count_;
            p_accum_ = static_cast<Real>(0.f);
            p_count_ = static_cast<Real>(1.f);
            Real dp = p - p_target_;
            d0_ = d0_ + 0.5f*dp*delta_d_*d0_;
            return d0_;
        }

        std::uniform_real_distribution<Real> new_displ_distr(){
            Real linear_displ = new_displacement_magnitude();
            return std::uniform_real_distribution<Real>(-linear_displ, linear_displ);
        }

        void probability_aggregator(std::optional<Real> e_diff, Real kBt){
            p_accum_ += prob(e_diff, kBt);
            p_count_ += static_cast<Real>(1.);
        }

        void reset_target_acceptance_probability(Real p_target){ p_target_ = p_target; }
        [[nodiscard]] Real target_acceptance_probability() const { return p_target_; }


    };

    [[maybe_unused]] [[nodiscard]] static Real node_curvature_squared(Index node_id, Nodes const& nodes){
        //! `unit_bending_energy` corresponds to the [Helfrich bending energy](https://en.wikipedia.org/wiki/Elasticity_of_cell_membranes) with bending rigidity 1 and gaussian bending stiffness 0.
        /**
         * \f[
         *  \mathrm{unit\_bending\_energy} = \frac{1}{2} A_{\mathrm{node}} (2 H_{node})^2
         * \f]
         * where \f$ H_{node} \f$ is the mean curvature of the node given by:
         * \f[
         * H_{node}^2 = \frac{\vec{K}_{node}}{2A_{node}} \cdot \frac{\vec{K}_{node}}{2A_{node}}
         * \f],
         * with  \f$ \vec{K} \f$ denoting the Node::curvature_vector.
         * @see Node::curvature_vec Triangulation::update_bulk_node_geometry(Index)
         */
        const vec3<Real> K = nodes.curvature_vec(node_id);
        const Real C2  = K.norm_square();///(static_cast<Real>(4.f)*A*A);

        return C2;
    }

    [[maybe_unused]] [[nodiscard]] static Real node_unit_bending_energy(Index node_id, Nodes const& nodes){
        //! `unit_bending_energy` corresponds to the [Helfrich bending energy](https://en.wikipedia.org/wiki/Elasticity_of_cell_membranes) with bending rigidity 1 and gaussian bending stiffness 0.
        /**
         * \f[
         *  \mathrm{unit\_bending\_energy} = \frac{1}{2} A_{\mathrm{node}} (2 H_{node})^2
         * \f]
         * where \f$ H_{node} \f$ is the mean curvature of the node given by:
         * \f[
         * H_{node}^2 = \frac{\vec{K}_{node}}{2A_{node}} \cdot \frac{\vec{K}_{node}}{2A_{node}}
         * \f],
         * with  \f$ \vec{K} \f$ denoting the Node::curvature_vector.
         * @see Node::curvature_vec Triangulation::update_bulk_node_geometry(Index)
         */
        const vec3<Real> K = nodes.curvature_vec(node_id);
        const Real A = nodes.area(node_id);
        const Real e  = static_cast<Real>(0.5)*A*K.norm_square();

        return e;
    }

}


#endif
