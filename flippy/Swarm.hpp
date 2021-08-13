#ifndef SYNC_BACK_WATCHER_PY_SWARM_HPP
#define SYNC_BACK_WATCHER_PY_SWARM_HPP

#include <random>
#include "vec3.hpp"

namespace fp {

template<typename Real>
class Particle
{
public:
    Particle() = default;
    Particle(Real R_inp, vec3<Real> pos_inp)
            :R_(R_inp)
    {
        set_pos(pos_inp);
    }
    vec3<Real> pos() const { return pos_; }
    Real R() const { return R_; }
//    Real overlap_energy() const { return overlap_energy_; }
    void set_pos(vec3<Real> const& new_pos) { pos_ = new_pos; }
    void update_pos(vec3<Real> const& pos_change) { set_pos(pos_ + pos_change); }
//    void update_overlap_energy(Real const& new_overlap_energy) { overlap_energy_=new_overlap_energy; }
private:
    Real R_;
//    Real overlap_energy_;
    vec3<Real> pos_{};
};
template<typename Type, typename Index, int size>
class Swarm
{
public:
    std::array<Particle<Type>, size> data;

    Swarm() = default;

    Swarm(Type const& R_inp, vec3<Type> const& origin, Type const& inner_shell_distance, Type const& middle_phi)
            :
            R_(R_inp)
    {
        Type sin_th;
        Type th = M_PI/2;
        std::array<Type, size> phis{-M_PI, middle_phi, 0};
        Type distance = inner_shell_distance + R_;
        Index i = 0;
        for (auto& particle: data) {
            sin_th = sin(th);
            vec3<Type> r{cos(phis[i])*sin_th, sin(phis[i])*sin_th, cos(th)};
            particle = Particle<Type>(R_, origin + distance*r);
            ++i;
        }
    }

    Swarm(Type const& R_inp, std::array<vec3<Type>, size> const& positions)
            :R_(R_inp)
    {
        for (std::size_t i = 0; i<size; ++i) {
            data[i] = Particle<Type>(R_, positions[i]);
        }
    }
    Type R() const { return R_; }
    vec3<Type> pos(Index particle_id) const { return data[particle_id].pos(); }
    void set_pos(Index particle_id, vec3<Type> const& new_pos) { data[particle_id].set_pos(new_pos); }
    void update_pos(Index particle_id, vec3<Type> const& pos_change) { data[particle_id].update_pos(pos_change); }
    void update_overlap_energy(Index particle_id, Type const& new_overlap_energy)
    {
        data[particle_id].update_overlap_energy(new_overlap_energy);
    }

private:
    Type R_;
//    std::uniform_real_distribution<Type> uniform_01_distr;
//    std::mt19937_64 random_engine;

};
}
#endif //SYNC_BACK_WATCHER_PY_SWARM_HPP
