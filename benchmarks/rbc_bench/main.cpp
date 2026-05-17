#if ALLOCATIONTRACKING

#include <cstdint>
#include <cstdlib>
#include <iostream>

namespace {
struct AllocMetrics {
    std::int64_t  allocation     = 0;
    std::uint64_t allocated_size = 0;
    std::int64_t  de_allocation  = 0;
    std::int64_t  leaks          = 0;
};
} // namespace

static AllocMetrics AM;

void *operator new(size_t size) {
    AM.allocation     += 1;
    AM.allocated_size += size;
    //    std::cout<<"allocating "<<(long double)size/1000.<<" kB\n";
    return std::malloc(size);
}
void operator delete(void *memory) noexcept {
    AM.de_allocation += 1;
    std::free(memory);
}

void operator delete(void *memory, std::size_t) noexcept {
    AM.de_allocation += 1;
    std::free(memory);
}

#endif

#include <FlippyBenchmarkLogger.h>
#include <flippy.hpp>

namespace {
using fp::operator""_r;
struct EnergyParameters {
    fp::Real kappa, K_V, K_A, V_t, A_t;
};
using Triangulation = fp::Triangulation;
using MCU           = fp::MonteCarloUpdater<EnergyParameters, std::mt19937>;

// This is the energy function that is used by flippy's built-in updater to decide if a move was energetically favorable
// or not
fp::Real surface_energy(const fp::Node               &node,
                        const fp::Triangulation      &trg,
                        const EnergyParameters       &prms,
                        const std::vector<fp::Index> &changed_neighborhood) {
    fp::Real V  = trg.global_geometry().volume;
    fp::Real A  = trg.global_geometry().area;
    fp::Real dV = V - prms.V_t;
    fp::Real dA = A - prms.A_t;

    fp::Real E_bending = fp::changed_neighborhood_bending_energy(node, trg.nodes(), changed_neighborhood);

    fp::Real energy = prms.kappa * E_bending + prms.K_V * dV * dV / prms.V_t + prms.K_A * dA * dA / prms.A_t;
    return energy;
}
} // namespace

int main(int /*argc*/, char **argv) {
    char LOG_FILE[300]{};
    {
        // triangulation iteration number of nodes
        auto     json_path = std::string(argv[1]);
        fp::Json config    = fp::json_read(json_path);
        // N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_trng
        unsigned int n_triang     = config["n_triang"].get<unsigned int>();
        fp::Real     l_min        = config["l_min"].get<fp::Real>();
        fp::Real     kappa        = config["kappa"].get<fp::Real>(); /*kBT*/
        fp::Real     K_A          = config["K_A"].get<fp::Real>();   /*kBT/volume*/
        fp::Real     K_V          = config["K_V"].get<fp::Real>();   /*kBT/area*/
        fp::Real     red_vol      = config["red_vol"].get<fp::Real>();
        int          max_mc_steps = config["max_mc_steps"].get<int>();
        std::string  save_dir     = config["save_dir"];

        // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial
        // bond length is close to minimal. This formula is derived from the equidistant sub-triangulation of an
        // icosahedron, where geodesic distances are used as a distance measure.
        fp::Real R     = fp::min_radius_with_non_overlapping_beads(l_min, n_triang);
        fp::Real l_max = 1.8 * l_min; // if you make l_max closer to l_min
        // bond_flip acceptance rate will go down
        fp::Real r_Verlet = 2 * l_max;

        fp::Real         Vt = red_vol * fp::sphere_vol(R);
        fp::Real         At = fp::sphere_area(R);
        EnergyParameters prms{.kappa = kappa, .K_V = K_V, .K_A = K_A, .V_t = Vt, .A_t = At};
        // side length of a voxel from which the displacement of the node is drawn
        fp::Real linear_displ = l_min / 8;
        // max number of iteration steps (depending on the strength of your CPU, this should take anywhere from a couple
        // of seconds to a couple of minutes

        std::random_device random_number_generator_seed;
        // create a random number generator and seed it with the current time
        std::mt19937 rng(random_number_generator_seed());

        auto guv        = Triangulation::make_spherical_triangulation(n_triang, R, r_Verlet);
        auto mc_updater = MCU(guv, prms, surface_energy, rng, l_min, l_max);

        fp::vec3<fp::Real>                       displ{};
        std::uniform_real_distribution<fp::Real> displ_distr(-linear_displ, linear_displ);

        // squish the sphere in the z-direction to break the initial symmetry. This speeds up the convergence to a
        // biconcave shape greatly.
        guv.scale_node_coordinates(1, 1, 0.8);

        fp::Json data_init = guv.make_egg_data();
        fp::json_dump(save_dir + "test_run_init", data_init);

        auto globals_saver = [&](const Triangulation &trg) {
            std::unordered_map<std::string, std::string> out;

            out["volume"]     = std::to_string(trg.global_geometry().volume);
            out["area"]       = std::to_string(trg.global_geometry().area);
            out["volume_rel"] = std::to_string(trg.global_geometry().volume / Vt);
            out["area_rel"]   = std::to_string(trg.global_geometry().area / At);
            return out;
        };

        auto stream_particle = [&](const fp::Node &node, const fp::Triangulation &) {
            std::string              s{"1"};
            fp::vec3<fp::Real>       pos = node.pos;
            static const std::string empty{" "};
            s += empty + std::to_string(pos.x) + empty + std::to_string(pos.y) + empty + std::to_string(pos.z) + empty +
                 std::to_string(node.curvature_vec.norm()) + "\n";
            return s;
        };
        std::vector<fp::experimental::xyzProperty> properties{
            fp::experimental::xyzProperty{
                .name = "species", .xyz_type = fp::experimental::XYZ_TYPE::S, .column_count = 1},
            fp::experimental::xyzProperty{.name = "pos", .xyz_type = fp::experimental::XYZ_TYPE::R, .column_count = 3},
            fp::experimental::xyzProperty{
                .name = "curv", .xyz_type = fp::experimental::XYZ_TYPE::R, .column_count = 1}};
        auto xyz_saver = fp::experimental::xyzDataSaver(save_dir + "data.xyz", properties, stream_particle);
        std::vector<fp::Index> shuffled_ids;
        shuffled_ids.reserve(guv.size());

        // create a vector that contains all node ids. We can shuffle this vector in each MC step to iterate randomly
        // through the nodes
        for (const auto &node : guv.nodes()) {
            shuffled_ids.push_back(node.id);
        }

        auto                    probability_target = 0.5_r;
        std::optional<fp::Real> e_diff             = 0_r;

        fp::Real V0            = guv.global_geometry().volume;
        fp::Real A0            = guv.global_geometry().area;
        auto     displ_updater = fp::DynamicDisplacementUpdater(linear_displ, probability_target);
        auto     logger        = FlippyBenchmarkLogger<MCU>::start_benchmark(mc_updater, guv);

        for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
            auto t      = static_cast<fp::Real>(mc_step);
            auto t_max  = static_cast<fp::Real>(max_mc_steps);
            prms.V_t    = fp::linear_adaptation(t, 0_r, 0.3_r * t_max, V0, Vt);
            prms.A_t    = fp::linear_adaptation(t, 0_r, 0.3_r * t_max, A0, At);
            fp::Real kT = fp::linear_adaptation(t, 0.7_r * t_max, 0.9_r * t_max, 1_r, 0.1_r);
            mc_updater.reset_kBT(kT);
            displ_distr = displ_updater.new_displ_distr();
            for (fp::Index node_id : shuffled_ids) {
                displ  = {.x = displ_distr(rng), .y = displ_distr(rng), .z = displ_distr(rng)};
                e_diff = mc_updater.move_MC_updater(guv[node_id], displ);
                displ_updater.probability_aggregator(e_diff, mc_updater.kBT());
            }
            fp::Real pt = fp::linear_adaptation(t, 0.7_r * t_max, 0.9_r * t_max, probability_target, 0.2_r);
            displ_updater.reset_target_acceptance_probability(pt);
            // then we shuffle the bead_ids
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            // then we loop through all of
            for (auto node_id : shuffled_ids) {
                // them again and try to flip their bonds
                mc_updater.flip_MC_updater(guv[node_id]);
            }
            if (mc_step % 1000 == 0) { xyz_saver.save_current_state_in_memory(guv, globals_saver); }
            if (mc_step % 10000 == 0) {
                xyz_saver.save_current_state_in_memory(guv, globals_saver);
                xyz_saver.flush_data_to_file();
            }
        }

        xyz_saver.save_current_state_in_memory(guv, globals_saver);
        xyz_saver.flush_data_to_file();

        std::string log_file_path = logger.log(config["log_dir"], config);

        const char *log_file_char = log_file_path.c_str();
        for (size_t i = 0; i < log_file_path.size(); ++i) {
            LOG_FILE[i] = log_file_char[i];
        }
        cutils::print(LOG_FILE);

        fp::Json data_final = guv.make_egg_data();
        // ATTENTION!!! this file will be saved in the same folder as the executable
        fp::json_dump(save_dir + "test_run_final", data_final);
    }

#if ALLOCATIONTRACKING
    AM.leaks = AM.allocation - AM.de_allocation;
    std::printf("total allocations: %5ld | total de_allocations: %5ld | leaks      %5ld\n",
                AM.allocation,
                AM.de_allocation,
                AM.leaks);
    std::printf("allocated size:    %5lo\n", AM.allocated_size);

    std::ofstream log(LOG_FILE, std::ios_base::app);
    log << "total_allocations: " << AM.allocation << '\n'
        << "total_de_allocations: " << AM.de_allocation << '\n'
        << "leaks: " << AM.leaks << '\n'
        << "allocated_size: " << AM.allocated_size << '\n'
        << "COMMENTS: |\n \n";
#endif
    return 0;
}
