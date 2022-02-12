#include <random> // needed for random displacement generation
#include <vector>
#include "flippy.hpp"


void append_XYZ_stream(std::stringstream &xyzStream, std::vector<int> const& node_ids, fp::Triangulation<double, int> const& guv ){

    xyzStream<<(node_ids.size())<<'\n';
    xyzStream<<"Atoms\n";

    for(int node_id: node_ids){
        if(node_id==0){
            xyzStream<<"0 "<<guv[node_id].pos.x<<' '<<guv[node_id].pos.y<<' '<<guv[node_id].pos.z<<'\n';
        }else{
            xyzStream<<"1 "<<guv[node_id].pos.x<<' '<<guv[node_id].pos.y<<' '<<guv[node_id].pos.z<<'\n';
        }
    }

};

double sphere_vol(double R){return 4./3. * M_PI *R*R*R;}
double sphere_area(double R){return 4. * M_PI *R*R;}

struct SimulationParameters{ // a data structure that can hold all simulation parameter, s.t. they can be easily passed to functions
    double bending_rigidity, K_V, K_A, R, V_t, A_t, dx, l_max_square;
    int n_triang, t_max;
}tension;

template<typename Real, typename Index>
Real surface_energy_area_volume_ensemble(fp::Node<Real, Index> const& node, fp::Triangulation<Real, Index> const& triangulation, SimulationParameters const& prms){
    Real V = triangulation.global_geometry().volume;
    Real A = triangulation.global_geometry().area;
    Real dV = V-prms.V_t;
    Real dA = A-prms.A_t;
    Real e_tot = (prms.bending_rigidity/2.)*triangulation.global_geometry().dA_K2 + prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return e_tot;
}

int main(){
    fp::print("starting"); // write the string "starting" to the standard out. fp is flippy's namespace and print is a built-in print function
    fp::Timer timer; //setting up a timer that will print to console how long the simulation took
    int n_triang = 8; // triangulation iteration number of nodes N_node=12+30*n_triang+20*n_triang(n_triang-1)/2
    double l_min = 2;
    double R = 1.5/sin(asin(1./(2*sin(2.*M_PI/5.)))/(n_triang+1.)); // estimate of a typical bond length in the initial triangulation. This formula is derived from equidistant subtriangulation of an icosahedron, where geodesic distances are used as a distance measure.
    double l_max = 2.5*l_min;
    double V0 = sphere_vol(R);
    SimulationParameters prms{
        .bending_rigidity = 10 /*kBT_*/, .K_V = 100 /*kBT_*/, .K_A=1000, .R=R /*a.u.*/,
        .V_t=V0, .A_t=sphere_area(R),
        .dx=l_min/8., // side length of a voxel from which the displacement of the node is drawn
        .n_triang=n_triang,
        .t_max=300 // max number of iteration steps (depending on the strength of your cpu, this should take anywhere from a couple of seconds to a couple of minutes
    };

    std::random_device random_number_generator_seed;
    std::mt19937 rng(random_number_generator_seed()); // create a random number generator and seed it with current time

    // All the flippy magic is happening on the following two lines
    fp::Triangulation<double, int> tr(prms.n_triang, prms.R, 2*l_max);
    fp::MonteCarloUpdater<double, int, SimulationParameters, std::mt19937, fp::SPHERICAL_TRIANGULATION> mc_updater(tr, prms, surface_energy_area_volume_ensemble<double, int>, rng, l_min, l_max);

    fp::vec3<double> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for a later use as a random direction vector
    std::uniform_real_distribution<double> dx_distr(-prms.dx, prms.dx); //define a distribution from which the small displacements in x y and z directions will be drawn

    tr.scale_node_coordinates(1,1,0.8); // squish the sphere in z direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly
    fp::Json data_init = {{"nodes", tr.make_egg_data()}};
    fp::json_dump("test_run_init", data_init);  // ATTENTION!!! this file will be saved in the same folder as the executable

    std::vector<int> shuffled_ids;
    shuffled_ids.reserve(tr.size());
    for(auto const& node: tr.nodes()){
        shuffled_ids.push_back(node.id);
    }
    std::stringstream xyz;
    int update_per_mc_step_max=200;
    double progress=0;

    for(int t=1; t<prms.t_max+1;++t){
        if(t%10==0){
            fp::print("V/V_t=",tr.global_geometry().volume/prms.V_t, "A/A_t=", tr.global_geometry().area/prms.A_t);
            tr.make_global_geometry();
            tr.make_verlet_list();
        }
        progress = t/((double)prms.t_max);
        if(t<prms.t_max/2){
            //for the first half of the simulation we gradually decrease the target volume. If we do not do this the volume term will dominate the energy, and we will get weird shapes!
            // since the MonteCarloUpdater takes a reference to the SimulationParameters struct, we can just change the parameters and flippy will have access to the updated version.
            prms.V_t = V0*(1.-0.4*(2*progress));
        }

        else if(t<(2*prms.t_max)/3){
            mc_updater.reset_kBT(1. - (progress-0.5)); //here we cool down the temperature to force the monte carlo simulation to go to the minimum configuration faster
        }
        for(int update_per_mc_step=0; update_per_mc_step<update_per_mc_step_max;++update_per_mc_step) {
            for (int node_id: shuffled_ids) {
                displ = {dx_distr(rng), dx_distr(rng), dx_distr(rng)};
                mc_updater.move_MC_updater(tr[node_id], displ);
            }
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            for (int node_id: shuffled_ids) {
                mc_updater.flip_MC_updater(tr[node_id]);
            }
        }
        append_XYZ_stream(xyz, shuffled_ids, tr);
    }
    std::ofstream oxyz("data.xyz");
    oxyz << xyz.str();
    oxyz.close();

    // MonteCarloUpdater counts the number of accepted and rejected moves, and it distinguishes between if a rejection
    // occurred because of the energy or the bond length constraint. We can use this to print a simple statistics here.
    // This will help us decide if our displacement size is too large for example.
    fp::print("percentage of failed moves: ",(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count())/((double)mc_updater.move_attempt_count()));
    fp::print("percentage of failed flips: ",(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count())/((double)mc_updater.flip_attempt_count()));

    fp::Json data_final = {{"nodes", tr.make_egg_data()}};
    fp::json_dump("test_run_final", data_final);  // ATTENTION!!! this file will be saved in the same folder as the executable
    timer.stop(); // strictly speaking this is not necessary the timer would stop and print the time automatically when it gets deleted
    return 0;
}