#if ALLOCATIONTRACKING
#include <iostream>


struct AllocMetrics{
  long allocation=0;
  unsigned long long allocated_size=0;
  long de_allocation=0;
  long leaks=0;
};
static AllocMetrics AM;
void* operator new(size_t size)
{
    AM.allocation+=1;
    AM.allocated_size+=size;
//    std::cout<<"allocating "<<(long double)size/1000.<<" kB\n";
    return malloc(size);
}
void operator delete(void* memory) noexcept{
    AM.de_allocation+=1;
    free(memory);
}

void operator delete(void* memory, std::size_t) noexcept{
    AM.de_allocation+=1;
    free(memory);
}

#endif

#include <flippy.hpp>
#include <FlippyBenchmarkLogger.h>

using REAL = float;
using INDEX = unsigned int;

REAL operator"" _r(long double value) {
    return static_cast<REAL>(value);
}
REAL operator"" _r(unsigned long long value) {
    return static_cast<REAL>(value);
}

struct EnergyParameters{REAL kappa, K_V, K_A, V_t, A_t;};
using Triangulation = fp::Triangulation<REAL, INDEX, fp::SPHERICAL_TRIANGULATION>;
using MCU = fp::MonteCarloUpdater<REAL, INDEX, EnergyParameters, std::mt19937, fp::SPHERICAL_TRIANGULATION>;

// This is the energy function that is used by flippy's built-in updater to decide if a move was energetically favorable or not
REAL surface_energy([[maybe_unused]]fp::Node<REAL, INDEX> const& node,
                      fp::Triangulation<REAL, INDEX> const& trg,
                      EnergyParameters const& prms,
                      std::vector<INDEX>const&
                      ){
    REAL V = trg.global_geometry().volume;
    REAL A = trg.global_geometry().area;
    REAL dV = V-prms.V_t;
    REAL dA = A-prms.A_t;
    REAL energy = prms.kappa*trg.global_geometry().unit_bending_energy +
                    prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return energy;
}

int main(int /*argc*/, char** argv) {
    char LOG_FILE[300]{};
    {
        // triangulation iteration number of nodes
        auto json_path = std::string(argv[1]);
        fp::Json config = fp::json_read(json_path);
        // N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_trng
        unsigned int n_triang = config["n_triang"].get<unsigned int>();
        REAL l_min = config["l_min"].get<REAL>();
        REAL kappa = config["kappa"].get<REAL>(); /*kBT*/
        REAL K_A = config["K_A"].get<REAL>(); /*kBT/volume*/
        REAL K_V = config["K_V"].get<REAL>(); /*kBT/area*/
        REAL red_vol = config["red_vol"].get<REAL>();
        int max_mc_steps = config["max_mc_steps"].get<int>();
        std::string save_dir = config["save_dir"];

        // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial bond length is close to minimal. This formula is derived from the equidistant sub-triangulation of an icosahedron, where geodesic distances are used as a distance measure.
        REAL R = fp::min_radius_with_non_overlapping_beads(l_min, n_triang);
        REAL l_max = 1.8_r * l_min; // if you make l_max closer to l_min
        // bond_flip acceptance rate will go down
        REAL r_Verlet = 2_r * l_max;

        REAL Vt = red_vol * fp::sphere_vol(R);
        REAL At = fp::sphere_area(R);
        EnergyParameters prms{.kappa= kappa, .K_V=K_V, .K_A=K_A, .V_t=Vt, .A_t=At};
        // side length of a voxel from which the displacement of the node is drawn
        REAL linear_displ = l_min / 8_r;
        // max number of iteration steps (depending on the strength of your CPU, this should take anywhere from a couple of seconds to a couple of minutes

        std::random_device random_number_generator_seed;
        // create a random number generator and seed it with the current time
        std::mt19937 rng(random_number_generator_seed());

        auto guv = Triangulation(n_triang, R, r_Verlet);
        auto mc_updater = MCU(guv, prms, surface_energy, rng, l_min, l_max);

        fp::vec3<REAL> displ{};
        std::uniform_real_distribution<REAL> displ_distr(-linear_displ, linear_displ);

        // squish the sphere in the z-direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly.
        guv.scale_node_coordinates(1_r, 1_r, 0.8_r);

        fp::Json data_init = guv.make_egg_data();
        fp::json_dump(save_dir+"test_run_init", data_init);

        std::vector<INDEX> shuffled_ids;
        shuffled_ids.reserve(guv.size());

        //create a vector that contains all node ids. We can shuffle this vector in each MC step to iterate randomly through the nodes
        for (auto const &node: guv.nodes()) {
            shuffled_ids.push_back(node.id);
        }

        REAL probability_target = 0.5f;
        REAL adaptation_magnitude = 0.1f;
        std::optional<REAL> e_diff = 0;

        REAL V0 = guv.global_geometry().volume;
        REAL A0 = guv.global_geometry().area;
        auto displ_updater = fp::DynamicDisplacementUpdater<REAL, INDEX>(linear_displ, adaptation_magnitude, probability_target);
        auto logger = FlippyBenchmarkLogger<REAL, INDEX, MCU>::start_benchmark(mc_updater, guv);

        for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
            auto t = static_cast<REAL>(mc_step);
            auto t_max = static_cast<REAL>(max_mc_steps);
            prms.V_t = fp::linear_adaptation(t, 0.f, 0.3f*t_max, V0, Vt);
            prms.A_t = fp::linear_adaptation(t, 0.f, 0.3f*t_max, A0, At);
            REAL kT = fp::linear_adaptation(t, 0.7f*t_max, 0.9f*t_max, 1.f, 0.1f);
            mc_updater.reset_kBT(kT);
            displ_distr = displ_updater.new_displ_distr();
            for (INDEX node_id: shuffled_ids) {
                displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
                e_diff = mc_updater.move_MC_updater(guv[node_id], displ);
                displ_updater.probability_aggregator(e_diff, mc_updater.kBT());
            }
            REAL pt = fp::linear_adaptation(t, 0.7f*t_max, 0.9f*t_max, probability_target, 0.2f);
            displ_updater.reset_target_acceptance_probability(pt);

            // then we shuffle the bead_ids
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            // then we loop through all of
            for (auto node_id: shuffled_ids) {
                // them again and try to flip their bonds
                mc_updater.flip_MC_updater(guv[node_id]);
            }
        }
        std::string log_file_path = logger.log(config["log_dir"], config);

        const char* log_file_char = log_file_path.c_str();
        for(size_t i=0; i<log_file_path.size(); ++i)
        {LOG_FILE[i]=log_file_char[i];}
        cutils::print(LOG_FILE);

        fp::Json data_final = guv.make_egg_data();
        // ATTENTION!!! this file will be saved in the same folder as the executable
        fp::json_dump(save_dir+"test_run_final", data_final);
    }

#if ALLOCATIONTRACKING
    AM.leaks = AM.allocation-AM.de_allocation;
    printf("total allocations: %5ld | total de_allocations: %5ld | leaks      %5ld\n",
           AM.allocation, AM.de_allocation, AM.leaks);
    printf("allocated size:    %5llo\n", AM.allocated_size);

    std::ofstream log(LOG_FILE, std::ios_base::app);
    log<<"total_allocations: "<<AM.allocation<<'\n'
       <<"total_de_allocations: "<<AM.de_allocation<<'\n'
       <<"leaks: "<<AM.leaks<<'\n'
       <<"allocated_size: "<<AM.allocated_size<<'\n'
       <<"COMMENTS: |\n \n";
#endif
    return 0;
}