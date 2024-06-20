#ifndef FLIPPY_SIM_UTILS_H
#define FLIPPY_SIM_UTILS_H
#include <custom_concepts.hpp>
#include <random>
#include <optional>

namespace fp{
    template<fp::floating_point_number Real>
    Real PI = static_cast<Real>(M_PI);

    template<fp::floating_point_number Real> Real sphere_vol(Real R){return static_cast<Real>(4./3.)* PI<Real> *R*R*R;}
    template<fp::floating_point_number Real> Real sphere_area(Real R){return static_cast<Real>(4.) * PI<Real> *R*R;}

    template<fp::floating_point_number Real> Real linear_adaptation(Real x, Real x_init, Real x_fin,
            Real val_init, Real val_fin)
     {
        if (x <= x_init) { return val_init; }
        if (x > x_fin) { return val_fin; }
//        if (std::abs(x_fin-x_init)<0.001){return val_init;}
        return val_init + (val_fin-val_init)/(x_fin-x_init)*(x-x_init);
    }

    template<fp::floating_point_number Real, fp::indexing_number Index>
    Index steps_for_fixed_speed_adaptation(Real x_init, Real val_init, Real val_fin, Real delta_val)
    {
        return static_cast<Index>(std::ceil((val_fin-val_init)/delta_val + x_init));
    }

    template<fp::floating_point_number Real, fp::indexing_number Index>
    Real min_radius_with_non_overlapping_beads(Real min_allowed_distance_betwean_bead_centers, Index sub_triangulation_iteration_count){
        Real l_min = min_allowed_distance_betwean_bead_centers;
        Real two = static_cast<Real>(2);
        Real one = static_cast<Real>(1);
        Real five = static_cast<Real>(5);
        Real n = static_cast<Real>(sub_triangulation_iteration_count);

        return l_min / (two * std::sin(std::asin(one / (two * std::sin(two * PI<Real> / five ))) / (n + one)));
    }



    template<fp::floating_point_number Real, fp::indexing_number Index>
    class DynamicDisplacementUpdater{
        Real d0_;
        Real delta_d_;
        Real p_target_;
        Real p_accum_{static_cast<Real>(0.f)};
        Real p_count_{static_cast<Real>(1.f)};

        Real prob(std::optional<Real> e_diff, Real kBT){ return (e_diff.has_value()?std::min(std::exp(e_diff.value()/kBT),1.f):0.f); }

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

}


#endif
