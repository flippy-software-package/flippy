#ifndef FLIPPY_TRIANGULATION_HPP
#define FLIPPY_TRIANGULATION_HPP

#include<optional>
#include "Nodes.hpp"
#include "vec3.hpp"
#include "utilities/utils.hpp"
#include "Triangulator.hpp"

/**
 * flippy's namespace.
 *
 * All class methods in flippy use similar prefix based nomenclature
 *
 * prefix | description
 * :-----:| :--------
 * calculate_| Indicates that a calculation will happen when the method is called, which might be expensive.
 * [action]_| [action]_ could be move_ or flip_ or any other descriptor. This prefixes indicate a state change and are usually expensive.
 * [no prefix] | usually signifies functions that return a constant reference to a private member (some people use get_ prefix for this).
 * Example: In the triangulation class `mass_center()` returns a const reference to the private member `mass_center.`
 */
namespace fp {

static constexpr int BOND_DONATION_CUTOFF = 4; // a node needs to have more than the cutoff number of bonds to be allowed to donate one

/**
 * A helper struct; keeps track of bond flips.
 * A bond flip can be unsuccessful, e.g. if the requested two nodes that are donating an edge already have too few edges.
 * If the flipp does happen then `common_nn_0` and `common_nn_1` record the ids of nodes that receive new common bond.
 *
 *```txt
 *    common_nn_0
 *     /         \
 *   /            \
 * node -------- nn_id
 *  \             /
 *   \           /
 *   common_nn_1
 *```
 *
 * */
template<integer_number Index>
struct BondFlipData
{
  bool flipped = false;
  Index common_nn_0 = -1;
  Index common_nn_1 = -1;
};

/**
 * A helper struct;  makes addition and subtraction on a ring easier.
 * */
template<integer_number Index>
struct Neighbors
{
  Index j_m_1{-1};  //neighbor j+1
  Index j_p_1{-1};  //neighbor j-1

  static Index plus_one(Index j, Index ring_size) { return ((j<ring_size - 1) ? j + 1 : (Index) 0); }
  static Index minus_one(Index j, Index ring_size) { return ((j==((Index) 0)) ? ring_size - 1 : j - 1); }

};
/**
 * A helper struct (template) that is used by the triangulation to pass data around in one convenient package.
 */
template<floating_point_number Real, integer_number Index>
struct Geometry
{
  Real area;
  Real volume;
  Real unit_bending_energy; //local area element times the square of the total curvature
  Geometry()
          :area(0.), volume(0.), unit_bending_energy(0.) { }
  explicit Geometry(Node<Real, Index> const& node)
          :area(node.area), volume(node.volume), unit_bending_energy(node.unit_bending_energy) { }
  Geometry(Real area_inp, Real volume_inp, Real unit_bending_energy_inp)
          :area(area_inp), volume(volume_inp), unit_bending_energy(unit_bending_energy_inp) { }
  friend Geometry<Real, Index> operator+(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
  {
      return Geometry<Real, Index>(lhs.area + rhs.area, lhs.volume + rhs.volume, lhs.unit_bending_energy + rhs.unit_bending_energy);
  }

  friend Geometry<Real, Index> operator-(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
  {
      return Geometry<Real, Index>(lhs.area - rhs.area, lhs.volume - rhs.volume, lhs.unit_bending_energy - rhs.unit_bending_energy);
  }

  void operator+=(Node<Real, Index> const& node)
  {
      area += node.area;
      volume += node.volume;
      unit_bending_energy += node.unit_bending_energy;
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

enum TriangulationType{
    SPHERICAL_TRIANGULATION, PLANAR_TRIANGULATION
};

/**
 * Implementation of Triangulation of closed two dimensional surfaces in 3D
 * See throughout the documentation the sources are referred to by numbers which can be looked up int the bibliography.
 *
 *
 * BIBLIOGRAPHY:
 *
 * [1] Guillaume Gueguen, Nicolas Destainville, and Manoel Manghi,
 * ‘Fluctuation Tension and Shape Transition of Vesicles: Renormalisation Calculations and Monte Carlo Simulations’,
 * Soft Matter, 13.36 (2017), 6100–6117 <https://doi.org/10.1039/C7SM01272A>.
 *
 * [2] Mark Meyer and others,
 * ‘Discrete Differential-Geometry Operators for Triangulated 2-Manifolds’,
 * in Visualization and Mathematics III, ed. by Hans-Christian Hege and Konrad Polthier,
 * Mathematics and Visualization (Berlin, Heidelberg: Springer, 2003),
 * pp. 35–57 <https://doi.org/10.1007/978-3-662-05105-4_2>.
 *
 */
template<floating_point_number Real, integer_number Index, TriangulationType triangulation_type=SPHERICAL_TRIANGULATION>
class Triangulation
{
private:
    explicit Triangulation(Real verlet_radius_inp)
    :global_geometry_(), verlet_radius(verlet_radius_inp){}
public:
    Triangulation() = default;
    //unit tested
    explicit Triangulation(Json const& nodes_input, Real verlet_radius_inp):Triangulation(verlet_radius_inp)
    {
        if constexpr(triangulation_type==SPHERICAL_TRIANGULATION) {
            nodes_ = Nodes<Real, Index>(nodes_input);
            all_nodes_are_bulk();
            initiate_advanced_geometry();
        }
        else{
            static_assert(triangulation_type==SPHERICAL_TRIANGULATION, "currently json initialization is only implemented for spherical triangulations!");
        }

    }

    Triangulation(Index n_nodes_iter, Real R_initial_input, Real verlet_radius_inp):Triangulation(verlet_radius_inp)
    {
        static_assert(triangulation_type==SPHERICAL_TRIANGULATION, "This initialization is intended for spherical triangulations");
        R_initial = R_initial_input;
        nodes_ = triangulate_sphere_nodes(n_nodes_iter);
        all_nodes_are_bulk();
        scale_all_nodes_to_R_init();
        orient_surface_of_a_sphere();
        initiate_advanced_geometry();
    }

    Triangulation(Index n_length, Index n_width, Real length, Real width, Real verlet_radius_inp):Triangulation(verlet_radius_inp)
    {
        static_assert(triangulation_type==PLANAR_TRIANGULATION, "This initialization is intended for planar triangulations");
        triangulate_planar_nodes(n_length, n_width, length, width);
        initiate_advanced_geometry();
    }

    void set_verlet_radius(Real R){
        verlet_radius = R;
        verlet_radius_squared = R*R;
    }

    //todo unittest
    void make_verlet_list()
    {
        for (auto& node: nodes_) {
            node.verlet_list.clear();
        }
        for (auto node_p = nodes_.begin(); node_p!=nodes_.end(); ++node_p) {
            for (auto other_node_p = nodes_.begin(); other_node_p!=node_p; ++other_node_p) {
                if ((node_p->pos - other_node_p->pos).norm_square()<verlet_radius_squared)
                {
                    node_p->verlet_list.push_back(other_node_p->id);
                    other_node_p->verlet_list.push_back(node_p->id);
                }
            }

        }
    }

    void translate_all_nodes(vec3<Real> const& translation_vector)
    {
        for (Index i = 0; i<nodes_.size(); ++i) { move_node(i, translation_vector); }
    }

    //unit tested
    vec3<Real> calculate_mass_center() const
    {
        vec3<Real> mass_center = vec3<Real>{0., 0., 0.};
        for (auto const& node : nodes_) { mass_center += node.pos; }
        mass_center = mass_center/nodes_.size();
        return mass_center;
    }

    //unit tested
    void move_node(Index node_id, vec3<Real> const& displacement_vector)
    {
        pre_update_geometry = get_two_ring_geometry(node_id);
        nodes_.displace(node_id, displacement_vector);
        update_two_ring_geometry(node_id);
        post_update_geometry = get_two_ring_geometry(node_id);
        update_global_geometry(pre_update_geometry, post_update_geometry);
    }

    // unit-tested
    void emplace_before(Index center_node_id, Index anchor_id, Index new_value)
    {
        /**
         * this method finds the anchor node in the nn_ids vector of the center_node
         * and uses Node classes own method emplace_nn_id to emplace the new_value
         * there (together with its distance to the center_node).
         * The body of this fuction looks like it does not guard against find returning
         * end() pointer, but this is taken care of in the emplace_nn_id method.
         */
        auto anchor_pos_ptr = std::find(nodes_[center_node_id].nn_ids.begin(),
                nodes_[center_node_id].nn_ids.end(), anchor_id);
        std::integral auto anchor_pos = (Index) (anchor_pos_ptr - nodes_[center_node_id].nn_ids.begin());
        nodes_[center_node_id].emplace_nn_id(new_value, nodes_[new_value].pos, anchor_pos);
    }

    //unit tested
    BondFlipData<Index> flip_bond(Index node_id, Index nn_id,
                                  Real min_bond_length_square,
                                  Real max_bond_length_square)
    {
        BondFlipData<Index> bfd{};
        if (nodes_.nn_ids(node_id).size()>BOND_DONATION_CUTOFF) {
            if (nodes_.nn_ids(nn_id).size()>BOND_DONATION_CUTOFF) {
                Neighbors<Index> common_nns = previous_and_next_neighbour_global_ids(node_id, nn_id);
                Real bond_length_square = (nodes_.pos(common_nns.j_m_1) - nodes_.pos(common_nns.j_p_1)).norm_square();
                if ((bond_length_square<max_bond_length_square) && (bond_length_square>min_bond_length_square)) {
                    if (common_neighbours(node_id, nn_id).size()==2) {
                    pre_update_geometry = calculate_diamond_geometry(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                    bfd = flip_bond_unchecked(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                        if (common_neighbours(bfd.common_nn_0, bfd.common_nn_1).size()==2) {
                    update_diamond_geometry(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                    post_update_geometry = calculate_diamond_geometry(node_id, nn_id, common_nns.j_m_1,
                            common_nns.j_p_1);
                    update_global_geometry(pre_update_geometry, post_update_geometry);
                        }
                        else {
                            flip_bond_unchecked(bfd.common_nn_0, bfd.common_nn_1, nn_id, node_id);
                            bfd.flipped = false;
                        }
                    }
                }
            }
        }
        return bfd;
    }

    //unit-tested
    void unflip_bond(Index node_id, Index nn_id, BondFlipData<Index> const& common_nns)
    {
        flip_bond_unchecked(common_nns.common_nn_0, common_nns.common_nn_1, nn_id, node_id);
//        copy_diamond_back(node_id, nn_id, cnn0_id, cnn1_id);
        update_diamond_geometry(node_id, nn_id, common_nns.common_nn_0, common_nns.common_nn_1);
        update_global_geometry(post_update_geometry, pre_update_geometry);
    }

    BondFlipData<Index> flip_bond_unchecked(Index node_id, Index nn_id,
                                            Index common_nn_j_m_1, Index common_nn_j_p_1)
    {
        emplace_before(common_nn_j_m_1, node_id, common_nn_j_p_1);
        emplace_before(common_nn_j_p_1, nn_id, common_nn_j_m_1);
        delete_connection_between_nodes_of_old_edge(node_id, nn_id);
        return {.flipped=true, .common_nn_0=common_nn_j_m_1, .common_nn_1=common_nn_j_p_1};
    }

    // unit-tested
    void update_bulk_node_geometry(Index node_id)
    {
        /**
         * calculates area volume and squared curvature integrated over the area, for the voronoi cell
         * associated to the node
         */
        update_nn_distance_vectors(node_id);

        Real area_sum = 0.;
        vec3<Real> face_normal_sum{0., 0., 0.}, local_curvature_vec{0., 0., 0.};
        vec3<Real> face_normal;
        std::integral auto nn_number = (Index) nodes_.nn_ids(node_id).size();
        Index j_p_1;

        Real face_area, face_normal_norm;
        vec3<Real> ljj_p_1, lij_p_1, lij;
        Real cot_at_j, cot_at_j_p_1;

        for (Index j = 0; j<nn_number; ++j) {
            //return j+1 element of ordered_nn_ids unless j has the last value then wrap around and return 0th element
            j_p_1 = Neighbors<Index>::plus_one(j,nn_number);

            lij = nodes_.nn_distances(node_id)[j];
            lij_p_1 = nodes_.nn_distances(node_id)[j_p_1];
            ljj_p_1 = lij_p_1 - lij;

            cot_at_j = cot_between_vectors(lij, (-1)*ljj_p_1);
            cot_at_j_p_1 = cot_between_vectors(lij_p_1, ljj_p_1);


            face_normal = lij.cross(lij_p_1); //nodes_.nn_distances(node_id)[j].cross(nodes_.nn_distances(node_id)[j_p_1]);
            face_normal_norm = face_normal.norm();
            face_area = mixed_area(lij, lij_p_1, 0.5*face_normal_norm, cot_at_j, cot_at_j_p_1);
            area_sum += face_area;
            face_normal_sum += face_area*face_normal/face_normal_norm;

            local_curvature_vec -= (cot_at_j_p_1*lij + cot_at_j*lij_p_1);
        }

        nodes_.set_area(node_id, area_sum);
        nodes_.set_volume(node_id, nodes_[node_id].pos.dot(face_normal_sum)/((Real) 3.)); // 18=3*6: 6 has the aforementioned justification. 3 is part of the formula for the tetrahedron volume
        nodes_.set_curvature_vec(node_id,  -local_curvature_vec/((Real) 2.*area_sum)); // 2 is part of the formula to calculate the local curvature I just did not divide the vector inside the loop
        nodes_.set_unit_bending_energy(node_id, local_curvature_vec.dot(local_curvature_vec)/((Real) 8.*area_sum)); // 8 is 2*4, where 4 is the square of the above two and the area in the denominator is what remains after canceling. 1/ comes from the pre-factor to bending energy

    };


    static std::tuple<Real, vec3<Real>> partial_voronoi_area_and_face_normal_of_node_in_a_triangle(vec3<Real> const& lij,
                                                                                                   vec3<Real> const& lij_p_1)
    {
        /** This function returns values of eqn.'s (82) & (84) from the paper [1]
         * Every node has its associated voronoi area and each voronoi area can be subdivided into parts that are
         * associated to each triangle that the node is part of. This function returns that sub-area and the face normal
         * of that triangle.
         */
        Real area, face_normal_norm;
        vec3<Real> un_noremd_face_normal;
        //precalculating this normal and its norm, will be needed in area calc. If all triangles are oriented as
        // right-handed, then this normal will point outwards
        un_noremd_face_normal = lij.cross(lij_p_1);
        face_normal_norm = un_noremd_face_normal.norm();
        area = mixed_area(lij, lij_p_1, face_normal_norm/2.);
        return std::make_tuple(area, un_noremd_face_normal);
    }

    static Real mixed_area(vec3<Real> const& lij, vec3<Real> const& lij_p_1, Real triangle_area, Real cot_at_j, Real cot_at_j_p_1){
        if ((cot_at_j>0.) && (cot_at_j_p_1>0.)) { // both angles at j and j+1 are smaller than 90 deg so the triangle can only be obtuse at the node
            if (lij.dot(lij_p_1)>0) { // cos at i is positive i.e. angle at i is not obtuse
                return (cot_at_j_p_1*lij.dot(lij) + cot_at_j*lij_p_1.dot(lij_p_1))/8.;
            }
            else {//obtuse at node i.
                return triangle_area/2.;
            }
        }
        else {//obtuse at node j or j+1.
            return triangle_area/4.;
        }

        }

//    unit tested
    [[deprecated("This function is deprecated and will be removed in a future release. mixed_area which does not take precalculated cotangents is performs expensive calculations use the use the alternative mixed_area function!")]]
    static Real mixed_area(vec3<Real> const& lij, vec3<Real> const& lij_p_1, Real const& triangle_area)
    {
        /** This function returns values of eqn.'s (82) (and two unnamed formulas in the following paragraph) from [1]
         *
         * Every node has its associated voronoi area and each voronoi area can be subdivided into parts that are
         * associated to each triangle that the node is part of. This function returns that sub-area. If the large triangle
         * (that the voronoi sub element is part of), is not obtuse. If it is obtuse, then the return value is either half,
         * or quarter of the large triangle. Depending where it is obtuse.
         *
         */
        vec3<Real> ljj_p_1 = lij_p_1 - lij;

        Real cot_at_j = cot_between_vectors(lij, (-1)*ljj_p_1);
        Real cot_at_j_p_1 = cot_between_vectors(lij_p_1, ljj_p_1);
        if ((cot_at_j>0.) && (cot_at_j_p_1>0.)) { // both angles at j and j+1 are smaller than 90 deg so the triangle can only be obtuse at the node
            if (lij.dot(lij_p_1)>0) { // cos at i is positive i.e. angle at i is not obtuse
                return (cot_at_j_p_1*lij.dot(lij) + cot_at_j*lij_p_1.dot(lij_p_1))/8.;
            }
            else {//obtuse at node i.
                return triangle_area/2.;
            }
        }
        else {//obtuse at node j or j+1.
            return triangle_area/4.;
        }

    }

    [[nodiscard]] Geometry<Real, Index> get_two_ring_geometry(Index node_id) const
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

    void update_two_ring_geometry(Index node_id)
    {
        /**
         * calculates area volume and squared curvature integrated over the area, for the  two-ring of the
         * associated to the node
         */
        update_bulk_node_geometry(node_id);
        for (auto nn_id: nodes_.nn_ids(node_id)) {
            update_bulk_node_geometry(nn_id);
        }
    };

    // unit-tested
    void scale_node_coordinates(Real x_stretch, Real y_stretch = 1, Real z_stretch = 1)
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
    [[nodiscard]] Geometry<Real, Index> calculate_diamond_geometry(Index node_id, Index nn_id,
                                                             Index cnn_0, Index cnn_1) const
    {
        /**
         * calculates area volume and squared curvature integrated over the area, for the diamond configuration of nodes
         * associated with a bond-flip from existing current_node_geometry
         */
        Geometry<Real, Index> diamond_geometry(nodes_[node_id]);
        diamond_geometry += nodes_[nn_id];
        diamond_geometry += nodes_[cnn_0];
        diamond_geometry += nodes_[cnn_1];
        return diamond_geometry;
    };

    //Todo unittest
    void update_diamond_geometry(Index node_id, Index nn_id, Index cnn_0, Index cnn_1)
    {
        /**
         * updates area volume and squared curvature integrated over the area, for the diamond configuration of nodes
         * associated with a bondflip
         */
        update_bulk_node_geometry(node_id);
        update_bulk_node_geometry(nn_id);
        update_bulk_node_geometry(cnn_0);
        update_bulk_node_geometry(cnn_1);
    };

    // Const Viewer Functions
    [[nodiscard]] Index size() const { return nodes_.size(); }
    const Node<Real, Index>& operator[](Index idx) const { return nodes_.data.at(idx); }
    const Nodes<Real, Index>& nodes() const { return nodes_; }
    [[nodiscard]] Json make_egg_data() const { return nodes_.make_data(); }
    [[nodiscard]] const Geometry<Real, Index>& global_geometry() const { return global_geometry_; }

    //Todo unittest
    void make_global_geometry()
    {
        const Geometry<Real, Index> empty{};
        global_geometry_ = empty;
//        for (auto const& node: nodes_.data) {
        for (auto node_id: bulk_nodes_ids) {
            update_bulk_node_geometry(node_id);
            update_global_geometry(empty, Geometry<Real, Index>(nodes_[node_id]));
        }
        for (auto node_id: boundary_nodes_ids) {
            update_boundary_node_geometry(node_id);
            update_global_geometry(empty, Geometry<Real, Index>(nodes_[node_id]));
        }
    }

    //Todo unittest
    void update_boundary_node_geometry(Index node_id){
        nodes_.set_area(node_id, 0.);
        nodes_.set_volume(node_id, 0.);
        nodes_.set_curvature_vec(node_id,  {0., 0., 0.});
        nodes_.set_unit_bending_energy(node_id, 0.);

    }

#ifdef TESTING_FLIPPY_TRIANGULATION_ndh6jclc0qnp274b
public:
#else
private:
#endif
    Real R_initial;
    Nodes<Real, Index> nodes_;
    std::vector<Index> boundary_nodes_ids;
    std::vector<Index> bulk_nodes_ids;
    Geometry<Real, Index> global_geometry_;
    Geometry<Real, Index> pre_update_geometry, post_update_geometry;
    mutable vec3<Real> l0_, l1_;
    Real verlet_radius{};
    Real verlet_radius_squared{};

    //unit tested
    void initiate_advanced_geometry(){
        initiate_distance_vectors();
        make_global_geometry();
        set_verlet_radius(verlet_radius);
        make_verlet_list();
    }

    //unit tested
    void scale_all_nodes_to_R_init()
    {
        static_assert(triangulation_type==SPHERICAL_TRIANGULATION, "This function is only well defined for a spherical triangulation");
        vec3<Real> diff;
        vec3<Real> mass_center = calculate_mass_center();
        for (Index i = 0; i<nodes_.size(); ++i) {
            diff = nodes_[i].pos - mass_center;
            diff.scale(R_initial/diff.norm());
            diff += mass_center;
            nodes_.set_pos(i, diff);
        }

    }

    void update_nn_distance_vectors(Index node_id)
    {
        /**
         *  This function calculates distance vectors from a node to all of it's neighbors. The directions of the
         *  distance vectors are (radiating outward) pointing from the node to neighbors.
         *  The function also preserves the order of neighbors. Meaning that the lilst of distance vectors is in the same
         *  order as th provided list of neighbor ids.
         */

//        nodes_[node_id].nn_distances.resize(nodes_[node_id].nn_ids.size());

        for (std::size_t i = 0; auto nn_id: nodes_.nn_ids(node_id)) {
            nodes_.set_nn_distance(node_id, i, nodes_.pos(nn_id) - nodes_.pos(node_id));
            ++i;
        }
    }

    Real cot_alphas_sum(Index node_id, Index nn_id, Index cnn_0, Index cnn_1) const
    {
        /**
         * given a node i and its neighbor j, they will share two common neighbor nodes p and m.
         * This function finds the angles at p & m opposite of i-j link.
         * This function implements the cot(alpha_ij) + cot(beta_ij) from fig. (6c) from [1].
         * The order of these neighbours does not matter for the correct sign of the angles.
         */

        l0_ = nodes_[node_id].pos - nodes_[cnn_0].pos;
        l1_ = nodes_[nn_id].pos - nodes_[cnn_0].pos;

        Real cot_sum = cot_between_vectors(l0_, l1_);
        l0_ = nodes_[node_id].pos - nodes_[cnn_1].pos;
        l1_ = nodes_[nn_id].pos - nodes_[cnn_1].pos;

        cot_sum += cot_between_vectors(l0_, l1_);
        return cot_sum;
    }

    static Real cot_between_vectors(vec3<Real> const& v1, vec3<Real> const& v2)
    {
        return v1.dot(v2)/(v1.cross(v2).norm());
    };

    //unit tested
    [[nodiscard]] std::vector<Index> order_nn_ids(Index node_id) const
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
        static_assert(triangulation_type==SPHERICAL_TRIANGULATION, "This function is only well defined for a spherical triangulation");
        std::vector<Index> nn_ids_temp;
        vec3<Real> li0, li1;
        vec3<Real> mass_center = calculate_mass_center();
        for (Index i = 0; i<nodes_.size(); ++i) { //ToDo modernize this loop
            nn_ids_temp = order_nn_ids(i);
            li0 = nodes_[nn_ids_temp[0]].pos - nodes_[i].pos;
            li1 = nodes_[nn_ids_temp[1]].pos - nodes_[i].pos;
            if ((li0.cross(li1)).dot(nodes_[i].pos - mass_center)<0) {
                std::reverse(nn_ids_temp.begin(), nn_ids_temp.end());
            }
            nodes_.set_nn_ids(i, nn_ids_temp);
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
    std::vector<Index> common_neighbours(Index node_id_0, Index node_id_1) const
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
    std::array<Index, 2> two_common_neighbours(Index node_id_0, Index node_id_1) const
    {
        std::array<Index, 2> res{-1, -1};
        //todo safe remove const& in the loop
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

    std::array<Index, 2> fast_two_common_neighbours(Index node_id_0, Index node_id_1) const
    {

        Index j = nodes_.find_nns_loc_idx(node_id_0, node_id_1);
        std::integral auto nn_number = (Index)nodes_.nn_ids(node_id_0).size();
        Index j_p_1 = Neighbors<Index>::plus_one(j, nn_number);
        Index j_m_1 = Neighbors<Index>::plus_one(j, nn_number);
        std::array<Index, 2> res{nodes_.nn_id(node_id_0,j_m_1),
                nodes_.nn_id(node_id_0,j_p_1)};
        return res;
    }

    std::array<Index, 2> two_common_neighbour_positions(Index node_id_0, Index node_id_1) const
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
    Neighbors<Index> previous_and_next_neighbour_local_ids(Index node_id, Index nn_id) const
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
    Neighbors<Index> previous_and_next_neighbour_global_ids(Index node_id, Index nn_id) const
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
    void delete_connection_between_nodes_of_old_edge(Index old_node_id0, Index old_node_id1)
    {
        nodes_[old_node_id0].pop_nn(old_node_id1);
        nodes_[old_node_id1].pop_nn(old_node_id0);
    }

    static Nodes<Real, Index> triangulate_sphere_nodes(Index n_iter){
        std::unordered_map<std::string,fp::implementation::SimpleNodeData<Real, Index>> simpleNodeData =
                fp::implementation::IcosahedronSubTriangulation<Real,Index>::make_corner_nodes();
        fp::implementation::IcosahedronSubTriangulation<Real,Index>::make_face_nodes(simpleNodeData, n_iter);

        Index nNewNodesOnEdge = n_iter - 1;
        Index nBulk = nNewNodesOnEdge*(nNewNodesOnEdge+1)/2;
        Index nNodes = fp::implementation::IcosahedronSubTriangulation<Real,Index>::N_ICOSA_NODEs
                + fp::implementation::IcosahedronSubTriangulation<Real,Index>::N_ICOSA_EDGEs*n_iter
                + fp::implementation::IcosahedronSubTriangulation<Real,Index>::N_ICOSA_FACEs*nBulk;
        std::vector<Node<Real, Index>> nodeData(nNodes);
        for(Index id; auto & nodeEl :simpleNodeData){
            id = nodeEl.second.id;
            nodeData[id].id = nodeEl.second.id;
            nodeData[id].pos = nodeEl.second.pos;
            for(auto const& hash: nodeEl.second.nn_hashes){
                nodeData[id].nn_ids.push_back(simpleNodeData[hash].id);
            }
        }
        return Nodes<Real, Index>(nodeData);
    }

    void triangulate_planar_nodes(Index n_length, Index n_width, Real length, Real width){
        Index N_nodes = n_length*n_width;
        fp::implementation::PlanarTriangulation<Real, Index> triang(n_length, n_width);
//        Nodes<Real, Index> bulk_nodes;
        Node<Real, Index> node;
        for(Index node_id=0; node_id<N_nodes; ++node_id){
            node.id = node_id;
            node.pos = fp::vec3<Real>{
                    triang.id_to_j(node_id)*length/n_length,
                    triang.id_to_i(node_id)*width/n_width,
                    0.
            };
            node.nn_ids = triang.nn_ids[node_id];
            nodes_.data.push_back(node);
            if(triang.is_bulk[node_id]){bulk_nodes_ids.push_back(node_id);}
            else{boundary_nodes_ids.push_back(node_id);}
        }
    }
    void all_nodes_are_bulk(){
        for(auto const& node: nodes_){
            bulk_nodes_ids.push_back(node.id);
        }
    }

};

}
#endif //FLIPPY_TRIANGULATION_HPP

