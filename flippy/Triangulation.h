/**
 * Implementation of Triangulation of closed two dimensional surfaces in 3D
 * See bibliography at the end of the file.
 */

#ifndef EXPERIMENT_H_TRIANGULATION_H_
#define EXPERIMENT_H_TRIANGULATION_H_

#include "external/delaunay-cpp-master/delaunay.hpp"
#include "external/json.hpp"

#include "Nodes.h"
#include "vec3.hpp"
#include "utilities/debug_utils.h"
#include "utilities/utils.h"

namespace fp {
using Json = nlohmann::json;

const int BOND_DONATION_CUTOFF = 5; // a node needs to have more than the cutoff number of bonds to be allowed to donate one
const int BOND_ACCEPTANCE_CUTOFF = 9; // a node needs to have less than the cutoff number of bonds to be allowed to accept one

template<typename Index>
struct BondFlipData
{
  bool flipped = false;
  Index common_nn_0 = -1;
  Index common_nn_1 = -1;
};

template<typename Index>
struct Neighbors
{
  Index j_m_1{-1};  //neighbor j+1
  Index j_p_1{-1};  //neighbor j-1

  static Index plus_one(Index j, Index ring_size) { return ((j<ring_size - 1) ? j + 1 : (Index) 0); }

  static Index minus_one(Index j, Index ring_size) { return ((j==((Index) 0)) ? ring_size - 1 : j - 1); }

};

template<typename Real, typename Index>
struct Geometry
{
  Real area;
  Real volume;
  Real dA_K2; //local area element times the square of the total curvature
  Geometry()
          :area(0.), volume(0.), dA_K2(0.) { }
  explicit Geometry(Node<Real, Index> const& node)
          :area(node.area), volume(node.volume), dA_K2(node.scaled_curvature_energy) { }
  Geometry(Real area_inp, Real volume_inp, Real dA_K2_inp)
          :area(area_inp), volume(volume_inp), dA_K2(dA_K2_inp) { }
  friend Geometry<Real, Index> operator+(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
  {
      return Geometry<Real, Index>(lhs.area + rhs.area, lhs.volume + rhs.volume, lhs.dA_K2 + rhs.dA_K2);
  }

  friend Geometry<Real, Index> operator-(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
  {
      return Geometry<Real, Index>(lhs.area - rhs.area, lhs.volume - rhs.volume, lhs.dA_K2 - rhs.dA_K2);
  }

  void operator+=(Node<Real, Index> const& node)
  {
      area += node.area;
      volume += node.volume;
      dA_K2 += node.scaled_curvature_energy;
  }

  friend void operator+=(Geometry<Real, Index>& lhs, Geometry<Real, Index> const& rhs)
  {
      lhs = lhs + rhs;
  }

  friend void operator-=(Geometry<Real, Index>& lhs, Geometry<Real, Index> const& rhs)
  {
      lhs = lhs - rhs;
  }
};

template<typename Real, typename Index>
class Triangulation
{
public:
    Triangulation() = default;
    //unit tested
    explicit Triangulation(Json const& nodes_input, Real const& verlet_radius)
            :nodes_(nodes_input, verlet_radius), mass_center_({0., 0., 0.}), global_geometry_()
    {
        initiate_simple_mass_center();
        initiate_distance_vectors();
        make_global_geometry();
        initiate_real_mass_center();
        make_verlet_list();
    }

    //unit tested
    Triangulation(Json const& nodes_input, Real R_initial_input, Real const& verlet_radius)
            :
            R_initial(R_initial_input), nodes_(nodes_input, verlet_radius), mass_center_({0., 0., 0.}),
            global_geometry_()
    {
        initiate_simple_mass_center();
        scale_all_nodes_to_R_init();
        orient_surface_of_a_sphere();
        initiate_distance_vectors(); //Todo if this is done before orient surface we can save time
        make_global_geometry();
        initiate_real_mass_center();
        make_verlet_list();
//        pre_move_nodes = nodes_;
    }

    Triangulation(Index const& nNodes, Real const& R_initial_input, Real const& verlet_radius_inp)
            :
            R_initial(R_initial_input), nodes_(), mass_center_({0., 0., 0.}), global_geometry_()
    {
        std::vector<vec3<Real>>
        nodePos = points_on_unit_sphere(nNodes);
        std::vector<std::vector<Index>> nnLists = find_nns(nodePos);
        std::vector<Node<Real, Index>>
        nodeData(nNodes);
        for (std::size_t idx = 0; Node<Real, Index>& node: nodeData) {
            node.pos = nodePos[idx];
            node.nn_ids = nnLists[idx];
            ++idx;
        }
        nodes_ = Nodes<Real, Index>(nodeData, verlet_radius_inp);

        initiate_simple_mass_center();
        scale_all_nodes_to_R_init();
        orient_surface_of_a_sphere();
        initiate_distance_vectors(); //Todo if this is done before orient surface we can save time
        make_global_geometry();
        initiate_real_mass_center();
        make_verlet_list();
    }

    vec3<Real> rS2(Real const& phi, Real const& theta)
    {
        return {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
    }

    template<typename T>
    int sgn(T val)
    {
        /*
         * taken from stack overflow:
         * https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
         */
        return (T(0)<val) - (val<T(0));
    }

    std::vector<vec3<Real>> points_on_unit_sphere(
    Index nNodes
    ){
        std::vector<vec3<Real>>
        positions(nNodes);
        for (Real phi = -1. + 1./(nNodes - 1.); auto& pos: positions) {
            pos = rS2(phi*(0.1 + 1.2*nNodes), M_PI*(1 - sgn(phi)*(1 - sqrt(1 - abs(phi))))/2.);
            phi += (2. - 2./(nNodes - 1.))/(nNodes - 1.);
        }
        return positions;
    }

    std::vector<delaunay::Point<float>> stereographic_projecttion(std::vector<vec3<Real>>
    const& nodePos){
        std::vector<delaunay::Point<float>> projectedPoints(nodePos.size());
        for (std::size_t idx = 0; auto& point: projectedPoints) {
            point.x = nodePos[idx][0]/(1 - nodePos[idx][2]);
            point.y = nodePos[idx][1]/(1. - nodePos[idx][2]);
            point.id = idx;
            ++idx;
        }
        return projectedPoints;
    }

    std::vector<std::vector<Index>> find_nns(std::vector<vec3<Real>>
    const& nodePos){
        std::vector<std::vector<Index>> nnList(nodePos.size());
        std::vector<delaunay::Point<float>> projectedPoints = stereographic_projecttion(nodePos);

        const auto triangulation = delaunay::triangulate(projectedPoints);
        std::vector<delaunay::Point<float>> plugProjectedPoints(5);
        std::copy(projectedPoints.end() - 5, projectedPoints.end(), plugProjectedPoints.begin());
        const auto plugTriangulation = delaunay::triangulate(projectedPoints);

        for (auto const& e : triangulation.edges) {
            print(e.p0.id, e.p1.id);
            nnList[e.p0.id].push_back(e.p1.id);
            nnList[e.p1.id].push_back(e.p0.id);
        }
        return nnList;
    }

    void make_verlet_list()
    {
        nodes_.make_verlet_list();
    }

    void translate_all_nodes(vec3<Real> const& translation_vector)
    {
        for (Index i = 0; i<nodes_.size(); ++i) { move_node(i, translation_vector); }
    }

    void initiate_real_mass_center()
    {
        mass_center_ = vec3<Real>{0, 0, 0};
        for (Index i = 0; i<nodes_.size(); ++i) { mass_center_ += nodes_[i].pos*nodes_[i].area; }
        mass_center_ = mass_center_/global_geometry_.area;
    }

    //unit tested
    void initiate_simple_mass_center()
    {
        mass_center_ = vec3<Real>{0, 0, 0};
        for (Index i = 0; i<nodes_.size(); ++i) { mass_center_ += nodes_[i].pos; }
        mass_center_ = mass_center_/nodes_.size();
    }

    //unit tested
    void move_node(Index node_id, vec3<Real> const& displacement_vector)
    {
//        copy_node_and_its_neighbours(node_id);
        pre_update_geometry = get_two_ring_geometry(node_id);
//        old_mass_center_=mass_center_;
        mass_center_ -= nodes_[node_id].pos*nodes_[node_id].area/global_geometry_.area;
        nodes_[node_id].pos += displacement_vector;
//        mass_center_ += displacement_vector/nodes_.size();
        update_two_ring_geometry(node_id);
        post_update_geometry = get_two_ring_geometry(node_id);
        update_global_geometry(pre_update_geometry, post_update_geometry);
        mass_center_ += nodes_[node_id].pos*nodes_[node_id].area
                /global_geometry_.area; //Todo make sure mass center is not needed in any geometry calculations

    }

    // Todo unittest
    void emplace_before(Index const& center_node_id, Index const& anchor_id, Index const& new_value)
    {
        auto anchor_pos_ptr = std::find(nodes_[center_node_id].nn_ids.begin(), nodes_[center_node_id].nn_ids.end(),
                anchor_id);
        auto anchor_pos = (Index) (anchor_pos_ptr - nodes_[center_node_id].nn_ids.begin());
        nodes_[center_node_id].emplace_nn_id(new_value, nodes_[new_value].pos, anchor_pos);
    }

    //Todo unittest
    BondFlipData<Index> flip_bond(Index const& node_id, Index const& nn_id, Real const& min_bond_length_square,
                                  Real const& max_bond_length_square)
    {
        BondFlipData<Index> bfd{};
        if (nodes_[node_id].nn_ids.size()>BOND_DONATION_CUTOFF) {
            if (nodes_[nn_id].nn_ids.size()>BOND_DONATION_CUTOFF) {
                Neighbors<Index> common_nns = previous_and_next_neighbour_global_ids(node_id, nn_id);
                Real bond_length_square = (nodes_[common_nns.j_m_1].pos - nodes_[common_nns.j_p_1].pos).norm_square();
                if ((bond_length_square<max_bond_length_square) && (bond_length_square>min_bond_length_square)) {
//                    if (common_neighbours(node_id, nn_id).size()==2) {

                    pre_update_geometry = get_diamond_geometry(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                    bfd = make_the_flip(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
//                        if (common_neighbours(bfd.common_nn_0, bfd.common_nn_1).size()==2) {
                    update_diamond_geometry(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                    post_update_geometry = get_diamond_geometry(node_id, nn_id, common_nns.j_m_1,
                            common_nns.j_p_1);
                    update_global_geometry(pre_update_geometry, post_update_geometry);
//                        }
//                        else {
//                            make_the_flip(bfd.common_nn_0, bfd.common_nn_1, nn_id, node_id);
//                            bfd.flipped = false;
//                        }
//                    }
                }
            }
        }
        return bfd;
    }

//    void unflip_bond(Index const& node_id, Index const& nn_id, Index const& cnn0_id, Index const& cnn1_id)
    void unflip_bond(Index const& node_id, Index const& nn_id, BondFlipData<Index> const& common_nns)
    {
        make_the_flip(common_nns.common_nn_0, common_nns.common_nn_1, nn_id, node_id);
//        copy_diamond_back(node_id, nn_id, cnn0_id, cnn1_id);
        update_diamond_geometry(node_id, nn_id, common_nns.common_nn_0, common_nns.common_nn_1);
        update_global_geometry(post_update_geometry, pre_update_geometry);
    }

    BondFlipData<Index> make_the_flip(Index const& node_id, Index const& nn_id,
                                      Index const& common_nn_j_m_1, Index const& common_nn_j_p_1)
    {
        emplace_before(common_nn_j_m_1, node_id, common_nn_j_p_1);
        emplace_before(common_nn_j_p_1, nn_id, common_nn_j_m_1);
        delete_connection_between_nodes_of_old_edge(node_id, nn_id);
        return {.flipped=true, .common_nn_0=common_nn_j_m_1, .common_nn_1=common_nn_j_p_1};
    }

    // unit-tested
    void update_node_geometry(Index const& node_id)
    {
        /**
         * calculates area volume and squared curvature integrated over the area, for the voronoi cell
         * associated to the node
         */
//        std::vector<Index> const& ordered_nn_ids = nodes_[node_id].nn_ids;
        update_nn_distance_vectors(node_id);
//        nodes_[node_id].nn_distances = all_nn_distance_vectors(node_id);

        Real area_sum = 0.;
        vec3<Real> face_normal_sum{0., 0., 0.}, local_curvature_vec{0., 0., 0.};
        vec3<Real> face_normal;
        auto nn_number = (Index) nodes_[node_id].nn_ids.size();
        Index j_p_1;
//        Real area;
//		vec3<Real> noremd_face_normal;
        for (Index j = 0; j<nn_number; ++j) {
            //return j+1 element of ordered_nn_ids unless j has the last value then wrap around and return 0th element
            j_p_1 = ((j<nn_number - 1) ? j + 1 : 0);
            face_normal = nodes_[node_id].nn_distances[j].cross(nodes_[node_id].nn_distances[j_p_1]);
            area_sum += face_normal.norm();
            face_normal_sum += face_normal;
            local_curvature_vec = local_curvature_vec + cot_alphas_sum(node_id, nodes_[node_id].nn_ids[j])
                    *nodes_[node_id]
                            .nn_distances[j]; //todo (speed) this is still too slow cos alphas are being over-calculated
        }
        // in all following cases 6=2*3; 2 comes from dividing face normal norm by 2 to get the right area and 3 from distributing the area over nodes
        area_sum = area_sum/((Real) 6.);
        nodes_[node_id].area = area_sum;
        nodes_[node_id].volume = nodes_[node_id].pos.dot(face_normal_sum)
                /((Real) 18.); // 18=3*6: 6 has the aforementioned justification. 3 is part of the formula for the tetrahedron volume
//        nodes_[node_id].volume = nodes_[node_id].pos.dot(face_normal_sum)/((Real) 3.); // 18=3*6: 6 has the aforementioned justification. 3 is part of the formula for the tetrahedron volume
        nodes_[node_id].curvature_vec = local_curvature_vec/((Real) 2.
                *area_sum); // 2 is part of the formula to calculate the local curvature I just did not divide the vector inside the loop
        nodes_[node_id].scaled_curvature_energy = local_curvature_vec.dot(local_curvature_vec)/((Real) 4.
                *area_sum); // 4 is the square of the above two and the area in the denominator is what remains after canceling

    };

    [[nodiscard]] Geometry<Real, Index> calculate_two_ring_geometry(Index const& node_id) const
    {
        /**
        * calculates area volume and squared curvature integrated over the area, for the  two-ring of the
        * associated to the node, using the stored current_node_geometry
        */
        Geometry<Real, Index> trg = node_geometry(node_id);
        for (auto const& nn_id: nodes_[node_id].nn_ids) {
            trg += node_geometry(nn_id);
        }
        return trg;
    }

    [[nodiscard]] Geometry<Real, Index> get_two_ring_geometry(Index const& node_id) const
    {
        /**
        * calculates area volume and squared curvature integrated over the area, for the  two-ring of the
        * associated to the node, using the stored current_node_geometry
        */
        Geometry<Real, Index> trg(nodes_[node_id]);
        for (auto const& nn_id: nodes_[node_id].nn_ids) {
            trg += nodes_[nn_id];
        }
        return trg;
    }

    void update_two_ring_geometry(Index const& node_id)
    {
        /**
         * calculates area volume and squared curvature integrated over the area, for the  two-ring of the
         * associated to the node
         */
        update_node_geometry(node_id);
        for (auto const& nn_id: nodes_[node_id].nn_ids) {
            update_node_geometry(nn_id);
        }
    };

    // unit-tested
    void ellipse_fy_cell(Real const& x_stretch, Real const& y_stretch = 1, Real const& z_stretch = 1)
    {
        vec3<Real> displ = {0, 0, 0};
        for (auto& node: nodes_.data) {
            displ[0] = node.pos[0]*(x_stretch - 1);
            displ[1] = node.pos[1]*(y_stretch - 1);
            displ[2] = node.pos[2]*(z_stretch - 1);
            move_node(node.id, displ);
        }
    }

    //Todo unittest
    [[nodiscard]] Geometry<Real, Index> get_diamond_geometry(Index const& node_id, Index const& nn_id,
                                                             Index const& cnn_0, Index const& cnn_1) const
    {
        /**
         * calculates area volume and squared curvature integrated over the area, for the diamond configuration of nodes
         * associated with a bond-flip from existing current_node_geometry
         */
        Geometry<Real, Index> diamond_geometry(nodes_[node_id]);
        diamond_geometry += nodes_[nn_id];
        diamond_geometry += nodes_[cnn_0];
        diamond_geometry += nodes_[cnn_1];
//        diamond_geometry += Geometry<Real, Index>(nodes_[nn_id]);
//        diamond_geometry += Geometry<Real, Index>(nodes_[cnn_0]);
//        diamond_geometry += Geometry<Real, Index>(nodes_[cnn_1]);

        return diamond_geometry;
    };

    //Todo unittest
    void update_diamond_geometry(Index const& node_id, Index const& nn_id, Index const& cnn_0, Index const& cnn_1)
    {
        /**
         * updates area volume and squared curvature integrated over the area, for the diamond configuration of nodes
         * associated with a bondflip
         */
        update_node_geometry(node_id);
        update_node_geometry(nn_id);
        update_node_geometry(cnn_0);
        update_node_geometry(cnn_1);
    };

    // Const Viewer Functions
    [[nodiscard]] Index size() const { return nodes_.size(); }
    const Node<Real, Index>& operator[](Index const idx) const { return nodes_.data.at(idx); }
    const vec3<Real>& mass_center() const { return mass_center_; }
    const Nodes<Real, Index>& nodes() const { return nodes_; }
    [[nodiscard]] Json egg_data() const { return nodes_.make_data(); }
    [[nodiscard]] const Geometry<Real, Index>& global_geometry() const { return global_geometry_; }

private:
    Real R_initial;
    Nodes<Real, Index> nodes_;
    vec3<Real> mass_center_;
    Geometry<Real, Index> global_geometry_;
    Geometry<Real, Index> pre_update_geometry, post_update_geometry;
    mutable vec3<Real> l0_, l1_;

    //unit tested
    void scale_all_nodes_to_R_init()
    {
        for (Index i = 0; i<nodes_.size(); ++i) {
            nodes_[i].pos =
                    mass_center_ + R_initial*(nodes_[i].pos - mass_center_)/(nodes_[i].pos - mass_center_).norm();
        }
        initiate_simple_mass_center();
    }

    //unit tested
    std::vector<vec3<Real>> all_nn_distance_vectors(
    Index const& node_id
    ) const
    {
        /**
         *  This function calculates distance vectors from a node to all of it's neighbors. The directions of the
         *  distance vectors are (radiating outward) pointing from the node to neighbors.
         *  The function also preserves the order of neighbors. Meaning that the lilst of distance vectors is in the same
         *  order as th provided list of neighbor ids.
         */

        std::vector<vec3<Real>>
        res;
        res.reserve(nodes_[node_id].nn_ids.size());

        for (auto const& nn_id: nodes_[node_id].nn_ids) { res.push_back(nodes_[nn_id].pos - nodes_[node_id].pos); }
        return res;
    }

    void update_nn_distance_vectors(Index const& node_id)
    {
        /**
         *  This function calculates distance vectors from a node to all of it's neighbors. The directions of the
         *  distance vectors are (radiating outward) pointing from the node to neighbors.
         *  The function also preserves the order of neighbors. Meaning that the lilst of distance vectors is in the same
         *  order as th provided list of neighbor ids.
         */

//        nodes_[node_id].nn_distances.resize(nodes_[node_id].nn_ids.size());

        for (std::size_t i = 0; auto const& nn_id: nodes_[node_id].nn_ids) {
            nodes_[node_id].nn_distances[i] = (nodes_[nn_id].pos - nodes_[node_id].pos);
            ++i;
        }
    }

    //unit tested
    Real cot_alphas_sum(Index const& node_id, Index const& nn_id) const
    {
        /**
         * given a node i and its neighbor j, they will share two common neighbor nodes p and m.
         * This function finds the angles at p & m opposite of i-j link.
         * This function implements the cot(alpha_ij) + cot(beta_ij) from fig. (6c) from [1].
         * The order of these neighbours does not matter for the correct sign of the angles.
         */
        auto common_nn_ids = two_common_neighbours(node_id, nn_id);
        l0_ = -nodes_[node_id].get_distance_vector_to(common_nn_ids[0]);
        l1_ = nodes_[nn_id].pos - nodes_[common_nn_ids[0]].pos;
        Real cot_sum = cot_between_vectors(l0_, l1_);
        l0_ = -nodes_[node_id].get_distance_vector_to(common_nn_ids[1]);
        l1_ = nodes_[nn_id].pos - nodes_[common_nn_ids[1]].pos;
        cot_sum += cot_between_vectors(l0_, l1_);
        return cot_sum;
    }

    static Real cot_between_vectors(vec3<Real> const& v1, vec3<Real> const& v2)
    {
        return v1.dot(v2)/(v1.cross(v2).norm());
    };

    //unit tested
    [[nodiscard]] std::vector<Index> order_nn_ids(Index const& node_id) const
    {
        std::vector<Index> const& nn_ids = nodes_[node_id].nn_ids;
        auto common_nn_ids = two_common_neighbours(node_id, nn_ids[0]);
        std::vector<Index> ordered_nn_ids{common_nn_ids[0], nn_ids[0], common_nn_ids[1]};

        Index nn_id;
        for (Index i = 0; i<(Index) nodes_[node_id].nn_ids.size() - 3; ++i) {
            nn_id = ordered_nn_ids[ordered_nn_ids.size() - 1];
            common_nn_ids = two_common_neighbours(node_id, nn_id);
            if (is_member(ordered_nn_ids, common_nn_ids[0])) {
                ordered_nn_ids.push_back(common_nn_ids[1]);
            }
            else {
                ordered_nn_ids.push_back(common_nn_ids[0]);
            }
        }

        return ordered_nn_ids;
    }

    //unit tested
    void orient_surface_of_a_sphere()
    {
        /**
         * If the initial configuration is spherical, then this function can orient the surface, such
         * that all right handed cross products will point outwards.
         * And the nn_ids are ordered in a way that two successive nn_s j and j+1 will give a right handed
         * cross product. I.e l_i_j x l_i_jp1 points outwards.
         *
         * This operation is not idempotent in a strict sense, since it guarantees that the nn_ids are in
         * a correct cycle every time but not in the same strict order, they might differ by an even
         * permutation. I.e. the ordering {1,2,3,4,5,6} and {6,1,2,3,4,5} are equivalent results.
         */
        std::vector<Index> nn_ids_temp;
        vec3<Real> li0, li1;

        for (Index i = 0; i<nodes_.size(); ++i) {
            nn_ids_temp = order_nn_ids(i);
            li0 = nodes_[nn_ids_temp[0]].pos - nodes_[i].pos;
            li1 = nodes_[nn_ids_temp[1]].pos - nodes_[i].pos;
            if ((li0.cross(li1)).dot(nodes_[i].pos - mass_center_)<0) {
                std::reverse(nn_ids_temp.begin(), nn_ids_temp.end());
            }
            nodes_[i].nn_ids = nn_ids_temp;
        }
    }

    // Todo unittest
    void initiate_distance_vectors()
    {
        for (Node<Real, Index>& node: nodes_.data) {
            node.nn_distances.resize(node.nn_ids.size());
            update_nn_distance_vectors(node.id);
        }
    }

    //unit tested
    std::vector<Index> common_neighbours(Index const& node_id_0, Index const& node_id_1) const
    {
        std::vector<Index> res;
        res.reserve(2);
        std::vector<Index> nn_ids0 = nodes_[node_id_0].nn_ids;
        std::vector<Index> nn_ids1 = nodes_[node_id_1].nn_ids;
        std::sort(nn_ids0.begin(), nn_ids0.end());
        std::sort(nn_ids1.begin(), nn_ids1.end());
        std::set_intersection(nn_ids0.begin(), nn_ids0.end(),
                nn_ids1.begin(), nn_ids1.end(),
                std::back_inserter(res));
        return res;
    }

    //unit tested
    std::array<Index, 2> two_common_neighbours(Index const& node_id_0, Index const& node_id_1) const
    {
        std::array<Index, 2> res{-1, -1};

        for (auto res_p = res.begin(); auto const& n0_nn_id: nodes_[node_id_0].nn_ids) {
            if (res_p==res.end()) { break; }
            else {
                if (is_member(nodes_[node_id_1].nn_ids, n0_nn_id)) {
                    *res_p = n0_nn_id;
                    ++res_p;
                }
            }
        }
        return res;
    }

    std::array<Index, 2> two_common_neighbour_positions(Index const& node_id_0, Index const& node_id_1) const
    {
        std::array<Index, 2> res{-1, -1};
        short counter = 0;
        for (auto const& n0_nn_id: nodes_[node_id_0].nn_ids) {
            if (counter==2) { break; }
            else {
                auto pos = std::find(nodes_[node_id_1].nn_ids.begin(), nodes_[node_id_1].nn_ids.end(), n0_nn_id);
                if (pos!=nodes_[node_id_1].nn_ids.end()) {
                    res[counter] = (Index) (pos - nodes_[node_id_1].nn_ids.begin());
                    ++counter;
                }
            }
        }
        return res;
    }

    //Todo unittest
    //unit tested
    Neighbors<Index> previous_and_next_neighbour_local_ids(Index const& node_id, Index const& nn_id) const
    {
        /**
         *        j+1
         *      /   \
         *     i-----j
         *     \    /
         *     	j-1
         *     	given i and j this function finds the local ids of j-1 and j+1 nodes and returns them IN THAT ORDER;
         *     	This function relies on the fact that i & j are neighbours and will throw a nasty runtime error if they are
         *     	not
         */
        auto const& nn_ids_view = nodes_[node_id].nn_ids;
        auto const local_nn_id = (Index) (std::find(nn_ids_view.begin(), nn_ids_view.end(), nn_id)
                - nn_ids_view.begin());
        auto const nn_number = (Index) nn_ids_view.size();
        return {.j_m_1= Neighbors<Index>::minus_one(local_nn_id, nn_number),
                .j_p_1 = Neighbors<Index>::plus_one(local_nn_id, nn_number)};
    }

    //unit tested
    Neighbors<Index> previous_and_next_neighbour_global_ids(Index const& node_id, Index const& nn_id) const
    {
        /**
         *        j+1
         *      /   \
         *     i-----j
         *     \    /
         *     	j-1
         *     	given i and j this function finds the global ids of j-1 and j+1 nodes and returns them IN THAT ORDER;
         *     	This function relies on the fact that i & j are neighbours and will throw a nasty runtime error if they are
         *     	not
         */
        auto const& nn_ids_view = nodes_[node_id].nn_ids;
        Neighbors<Index> neighbors = previous_and_next_neighbour_local_ids(node_id, nn_id);
        return {.j_m_1=nn_ids_view[neighbors.j_m_1], .j_p_1=nn_ids_view[neighbors.j_p_1]};
    }

    void update_global_geometry(Geometry<Real, Index> const& lg_old, Geometry<Real, Index> const& lg_new)
    {
        global_geometry_ += lg_new - lg_old;
    }

    // Todo unittest
    void delete_connection_between_nodes_of_old_edge(int const& old_node_id0, int const& old_node_id1)
    {
        nodes_[old_node_id0].pop_nn(old_node_id1);
        nodes_[old_node_id1].pop_nn(old_node_id0);
    }

    //Todo unittest
    void make_global_geometry()
    {
        const Geometry<Real, Index> empty{};
        global_geometry_ = empty;
        for (auto const& node: nodes_.data) {
            update_node_geometry(node.id);
            update_global_geometry(empty, Geometry<Real, Index>(node));
        }
    }

};
}
#endif //EXPERIMENT_H_TRIANGULATION_H_

/**
 * BIBLIOGRAPHY:
 *
 * [1] Guillaume Gueguen, Nicolas Destainville, and Manoel Manghi,
 * ‘Fluctuation Tension and Shape Transition of Vesicles: Renormalisation Calculations and Monte Carlo Simulations’,
 * Soft Matter, 13.36 (2017), 6100–6117 <https://doi.org/10.1039/C7SM01272A>.
 * [2] Mark Meyer and others,
 * ‘Discrete Differential-Geometry Operators for Triangulated 2-Manifolds’,
 * in Visualization and Mathematics III, ed. by Hans-Christian Hege and Konrad Polthier,
 * Mathematics and Visualization (Berlin, Heidelberg: Springer, 2003),
 * pp. 35–57 <https://doi.org/10.1007/978-3-662-05105-4_2>.
 *
 */