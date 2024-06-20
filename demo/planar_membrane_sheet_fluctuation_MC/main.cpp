#include <random> // needed for random displacement generation
#include <vector> // need for std::vector
#include <iostream> // needed for std::cout
#include "flippy.hpp"

// exyzStream is a class that can be used to stream data to a file in extended xyz format
// (see https://docs.ovito.org/reference/file_formats/input/xyz.html)
// Then the file can be read by ovito to visualize the data.
template<fp::floating_point_number Real, fp::indexing_number Index>
class exyzStream
{
    std::string xyzStream;
    std::vector<Index>* shuffled_ids_p;
    fp::Triangulation<Real,Index,fp::EXPERIMENTAL_PLANAR_TRIANGULATION>* triangulation;
    std::string size_string{};
    std::string properties_string{"Properties=species:S:1:pos:R:3\n"};
public:
    exyzStream()=default;
    exyzStream(std::vector<Index>* shuffled_ids_p_inp ,
               fp::Triangulation<Real, Index, fp::EXPERIMENTAL_PLANAR_TRIANGULATION>* guv_p_inp)
            : shuffled_ids_p(shuffled_ids_p_inp),
              triangulation(guv_p_inp),
              size_string(std::to_string(shuffled_ids_p->size())+'\n'){}


    void append_XYZ_stream_with_test_node(Index test_node_id){
        xyzStream.append(size_string);
        xyzStream.append(properties_string);
        auto test_node_neighbours = (*triangulation)[test_node_id].nn_ids;
        for(Index node_id: *shuffled_ids_p){
            if(node_id==test_node_id){
                xyzStream.append(stream_particle("11", (*triangulation)[node_id].pos));
            }
            else if(fp::is_member(test_node_neighbours, node_id)){
                xyzStream.append(stream_particle("12", (*triangulation)[node_id].pos));
            }
            else{
                xyzStream.append(stream_particle("1", (*triangulation)[node_id].pos));
            }
        }
    }
    void append_XYZ_stream(){
        xyzStream.append(size_string);
        xyzStream.append(properties_string);
        for(Index node_id: *shuffled_ids_p){
            xyzStream.append(stream_particle("1", (*triangulation)[node_id].pos));
        }
    }

    void streamXYZ(){
        std::ofstream data_file;
        data_file.open("data.xyz");
        data_file<<xyzStream;
    }

    static std::string stream_particle(std::string const& name, auto const& vec){
        std::string s;
        static const std::string empty{" "};
        s = name
            +empty+std::to_string(vec[0])
            +empty+std::to_string(vec[1])
            +empty+std::to_string(vec[2])
            +"\n";
        return s;
    }


};

struct EnergyParameters{double kappa, K_A, A_t;};

// This is the energy function that is used by flippy's built-in updater to decide if a move was energetically favorable or not
double surface_energy([[maybe_unused]]fp::Node<double, unsigned int> const& node,
                      fp::Triangulation<double, unsigned int, fp::EXPERIMENTAL_PLANAR_TRIANGULATION> const& trg,
                      EnergyParameters const& prms,
                      std::vector<unsigned int>&){
    double A = trg.global_geometry().area;
    double dA = A-prms.A_t;
    double energy = prms.kappa*trg.global_geometry().unit_bending_energy + prms.K_A*dA*dA/prms.A_t;
    return energy;
}

int main(){
    double l_min = 2;
    unsigned int n_x = 30;
    unsigned int n_y = 30;
    double non_overlap_stretch = 1.05 ; //This stretch factor makes sure that the nodes do not overlap and the bonds are long enouth, s.t. the triangles are not degenerate (i.e. all nodes are on the same line)
    double l_x = non_overlap_stretch *(n_x-1)*l_min;
    double l_y = non_overlap_stretch *(n_x-1)*l_min;
    double l_max = 2.*l_min; // if you make l_max closer to l_min bond_flip acceptance rate will go down. However if l_max is large enough that degenerate triangles can forme, then the simulation will give wrong results.
    double r_Verlet = 2*l_max;
    EnergyParameters prms{.kappa=2 /*kBT*/, .K_A=1000 /*kBT/volume*/, .A_t=l_x*l_y};
    double linear_displ=l_min/10.; // side length of a voxel from which the displacement of the node is drawn
    int max_mc_steps=2e5; // max number of iteration steps (depending on the strength of your CPU, this should take anywhere from a couple of seconds to a couple of minutes

    std::random_device random_number_generator_seed;
    auto seed = random_number_generator_seed();
    std::cout<<"Seed: "<<seed<<'\n';
    std::mt19937 rng(seed); // create a random number generator and seed it with the current time

    // All the flippy magic is happening on the following two lines
    fp::Triangulation<double, unsigned int, fp::EXPERIMENTAL_PLANAR_TRIANGULATION> planar_trg(n_x, n_y, l_x, l_y, r_Verlet);
    fp::MonteCarloUpdater<double, unsigned int, EnergyParameters, std::mt19937, fp::EXPERIMENTAL_PLANAR_TRIANGULATION> mc_updater(planar_trg, prms, surface_energy, rng, l_min, l_max);

    fp::vec3<double> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for later use as a random direction vector
    std::uniform_real_distribution<double> displ_distr(-linear_displ, linear_displ); //define a distribution from which the small displacements in x y and z directions will be drawn

    fp::Json data_init = planar_trg.make_egg_data();
    fp::json_dump("test_run_init", data_init);  // ATTENTION!!! this file will be saved in the same folder as the executable

    std::vector<unsigned int> shuffled_ids;
    shuffled_ids.reserve(planar_trg.size());
    for(auto const& node: planar_trg.nodes()){ shuffled_ids.push_back(node.id);} //create a vector that contains all node ids. We can shuffle this vector in each MC step to iterate randomly through the nodes
    exyzStream<double, unsigned int> xyzStream(&shuffled_ids, &planar_trg);
    xyzStream.append_XYZ_stream();
    xyzStream.streamXYZ();



    for(int mc_step=0; mc_step<max_mc_steps; ++mc_step){
        for (unsigned int node_id: shuffled_ids) { // we first loop through all the beads and move them
            displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
            mc_updater.move_MC_updater(planar_trg[node_id], displ);
        }

        std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng); // then we shuffle the bead_ids
        for (unsigned int node_id: shuffled_ids) { // then we loop through all of them again and try to flip their bonds
            mc_updater.flip_MC_updater(planar_trg[node_id]);
        }
        if(mc_step>=max_mc_steps/2){
            mc_updater.reset_kBT((1.-2.*(static_cast<double>(mc_step)/static_cast<double>(max_mc_steps)-0.5))); // this is a simple way to decrease the temperature of the system. We could also have used a more sophisticated annealing schedule, where we cycle the temperature.
        }
        if(mc_step%300==0){
            xyzStream.append_XYZ_stream();
            xyzStream.streamXYZ(); // ATTENTION!!! this file will be saved in the same folder as the executable
            std::cout<<"mc_step: "<<mc_step<<'\n';
            std::cout<<"Energy: "<<planar_trg.global_geometry().unit_bending_energy <<'\n';
            std::cout<<"-------------------------\n";
        }
    }
    xyzStream.streamXYZ();

    // MonteCarloUpdater counts the number of accepted and rejected moves, distinguishing whether a rejection occurred because of the energy or the bond length constraint.
    // We can use this to print simple statistics here. For example, this will help us decide if our displacement size is too large.
    std::cout<<"percentage of failed moves: "<<(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count())/((long double)mc_updater.move_attempt_count())<<'\n';
    std::cout<<"percentage of failed flips: "<<(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count())/((long double)mc_updater.flip_attempt_count())<<'\n';

    fp::Json data_final = planar_trg.make_egg_data();
    fp::json_dump("test_run_final", data_final);  // ATTENTION!!! this file will be saved in the same folder as the executable

    return 0;
}