#ifndef FLIPPY_TRIANGULATION_HPP
#define FLIPPY_TRIANGULATION_HPP
/**
 * @file
 * @brief This file contains the fp::Triangulation class and several related helper classes. This is the core of flippy.
 */
#include<optional>
#include <set>
#include "Nodes.hpp"
#include "vec3.hpp"
#include "utilities/utils.hpp"
#include "Triangulator.hpp"
#include "stlSerializer.hpp"

/**
 * @GlobalsStub
 * @{
 */
//! Redefinition of the `LONG_LONG_MAX` macro in case a specific compiler does not implement it. The largest number that fits in the `long long` type.
#ifndef LONG_LONG_MAX
#define LONG_LONG_MAX 9223372036854775807LL
#endif

/**
 * @brief Literal for a very large integral number.
 *
 * This number will be used for the default instantiation of variables/ data members that hold indices.
 * This is done to avoid instantiation with 0, which could also be a valid index.
 * This way, one can differentiate between unintentionally default instantiated value and a proper value of an index.
 * This will also mean that the use of such improperly instantiated indexer variables will cause an easier-to-identify error during runtime.
 * @see fp::BondFlipData
 * fp::Neighbors
 * fp::Triangulation::two_common_neighbours(Index, Index) const
 * fp::Triangulation::two_common_neighbour_positions(Index, Index) const
 */
#define VERY_LARGE_NUMBER_ LONG_LONG_MAX
/**@}*/

/**
 * flippy's namespace.
 */
namespace fp {
/**
 * @GlobalsStub
 * @{
 */
//! a node needs to have more than the cutoff number of bonds to be allowed to donate one
static constexpr int BOND_DONATION_CUTOFF = 4;
/**@}*/

//! A helper struct; keeps track of bond flips.
/**
 * A bond flip can be unsuccessful, e.g., if the requested two nodes that are donating an edge already have too few edges
 * (more details on bond flips and how they can fail are provided in the Triangulation.flip_bond(Index, Index, Real, Real) function).
 * If the flip does happen, then #flipped will be changed to true by the Triangulation.flip_bond(Index, Index, Real, Real) function and
 * #common_nn_0 and #common_nn_1 will record the ids of nodes that receive new common bond.
 * If the flip does not happen then the #common_nn_0 and #common_nn_1 data members will hold a #VERY_LARGE_NUMBER_ with which they were initiated.
 *
 * The initiation with #VERY_LARGE_NUMBER_ is done to avoid #common_nn_0 and #common_nn_1 being zero initiated.
 * This mechanism will lead to wrong behavior that is easier to identify during debugging in situations of un-careful usage.
 * Example of such un-careful usage is if the end-user does not check that the bonds were not flipped and tries to un-flip bonds.
 *
 *```txt
 *
 * State of the triangulation before and after the flip.
 *  before the flip                  after the flip
 *
 *      common nn 1                     common nn 1
 *      /          \                    /    |    \
 *     /            \                  /     |     \
 *   node --------- nn              node     |     nn
 *     \            /                 \      |     /
 *      \          /                   \     |    /
 *     common nn 0                     common nn 0
 *```
 *
 * @tparam Index @IndexStub
 * @see Triangulation.flip_bond(Index, Index, Real, Real)
 * @see Triangulation.unflip_bond(Index, Index, BondFlipData<Index> const&)
 * @see Triangulation.flip_bond_unchecked(Index, Index, Index, Index)
 */
template<indexing_number Index>
struct BondFlipData
{
  bool flipped = false; //!< track if the bond was flipped.
  Index common_nn_0 = static_cast<Index>(VERY_LARGE_NUMBER_); //!< Global id of a node that is supposed to receive a bond after a flip.
  Index common_nn_1 = static_cast<Index>(VERY_LARGE_NUMBER_); //!< Global id of a node that is supposed to receive a bond after a flip.
};

/**
 * A helper struct;  makes addition and subtraction on a ring easier.
 * Each fp::Node stores its next neighbors in a vector Node.nn_ids,
 * where the adjacent members in a vector are also next neighbors of each other.
 * This struct provides a safe way to always access the next or previous member of the nn_ids vector, even if a wraparound is necessary.
 * @see fp::Node.
 * @tparam Index @IndexStub
 */
template<indexing_number Index>
struct Neighbors
{
  //! neighbor j+1
  Index j_m_1{VERY_LARGE_NUMBER_};
  //! neighbor j-1
  Index j_p_1{VERY_LARGE_NUMBER_};

  /**
   *
   * @param j index of the j-th next neighbor
   * @param ring_size corresponds to the number of elements of Node.nn_ids
   * @return `j+1` if `j < ring_size - 1` and `0` otherwise. I.e., if `j` is the index of Node.nn_ids, its next neighbor will be stored in `j+1` element, unless `j` is the last element, then its next neighbor will be stored in the 0th element.
   */
  static Index plus_one(Index j, Index ring_size) { return ((j<ring_size - 1) ? j + 1 : (Index) 0); }
  /**
 *
 * @param j index of the j-th next neighbor
 * @param ring_size corresponds to the number of elements of Node.nn_ids
 * @return `j-1` if `j > 0` and `ring_size - 1` otherwise. I.e., if `j` is the index of Node.nn_ids, its previous next neighbor will be stored in `j-1` element, unless `j` is the 0th element, then its previous next neighbor will be stored in the last element.
 */
  static Index minus_one(Index j, Index ring_size) { return ((j==((Index) 0)) ? ring_size - 1 : j - 1); }
};
//! A helper struct. Used by the triangulation class to pass data around in one convenient package.
/**
 * Geometry is a struct that contains the usually needed geometric data in a triangulation.
 * This struct can hold such data for a single node, a collection of nodes, or the entire triangulation.
 * In the abstract, the fp::Geometry struct contains geometric data associated with some surface patch.
 * This is useful since we often need to aggregate information for a node and its neighboring nodes.
 *
 * The Data members of this struct are public, and thus,
 * it does not guarantee the correctness or consistency of the data it holds,
 * since it can be changed externally.
 * The geometry struct also provides overloaded arithmetic operators,
 * that are usually useful when aggregating geometric data from single nodes over a larger patch of the surface.
 * @tparam Real @RealStub
 * @tparam Index @IndexStub
 */
template<floating_point_number Real, indexing_number Index>
struct Geometry
{
  Real area; //!< Area of the patch. Sum over the associated areas of individual nodes that comprise the patch. (Compare to Node::area).
  Real volume; //!< Volume of the patch. Sum over the associated volumes of individual nodes that comprise the patch. (Compare to Node::volume).
  Real unit_bending_energy; //!< Volume of the patch. Sum over the associated unit bending energies of individual nodes that comprise the patch. (Compare to Node::unit_bending_energy).
  //! Default constructor, that zero initiates all the data members.
  Geometry() :area(0.), volume(0.), unit_bending_energy(0.) { }
  //! Construct from a node.
  explicit Geometry(Node<Real, Index> const& node)
  :area(node.area), volume(node.volume), unit_bending_energy(node.unit_bending_energy) {
  /**
   * @param node Initiates data members of the struct with the geometric values of this single node.
   */
  }
  //! Direct Constructor.
  Geometry(Real area_inp, Real volume_inp, Real unit_bending_energy_inp)
  :area(area_inp), volume(volume_inp), unit_bending_energy(unit_bending_energy_inp) {
  /**
   * Initiates data members from directly provided numeric values.
   * @param area_inp Input area.
   * @param volume_inp Input volume.
   * @param unit_bending_energy_inp Input unit bending energy.
   */
  }

  //! Overloaded addition operator.
  /**
   * Implements addition between two Geometries.
   * Each data member of *lhs* and *rhs* Geometries is added pairwise and
   * stored in the corresponding data member of the returned struct.
   * ```c++
   * //Example
   * Geometry<float, unsigned int> res = lhs + rhs;
   *
   * res.area == lhs.area + rhs.area; // true
   * res.volume == lhs.volume + rhs.volume; // true
   * res.unit_bending_energy == lhs.unit_bending_energy + rhs.unit_bending_energy; // true
   * ```
   * @param lhs
   * @param rhs
   * @return lhs+rhs
   */
  friend Geometry<Real, Index> operator+(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
  {

      return Geometry<Real, Index>(lhs.area + rhs.area, lhs.volume + rhs.volume, lhs.unit_bending_energy + rhs.unit_bending_energy);
  }
    //! Overloaded subtraction operator.
    /**
     * Implements subtraction between two Geometries.
     * Each data member of *lhs* and *rhs* Geometries is subtracted pairwise and
     * stored in the corresponding data member of the returned struct.
     * ```c++
     * //Example
     * Geometry<float, unsigned int> res = lhs - rhs;
     *
     * res.area == lhs.area - rhs.area; // true
     * res.volume == lhs.volume - rhs.volume; // true
     * res.unit_bending_energy == lhs.unit_bending_energy - rhs.unit_bending_energy; // true
     * ```
     * @param lhs
     * @param rhs
     * @return lhs-rhs
     */
  friend Geometry<Real, Index> operator-(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
  {
      return Geometry<Real, Index>(lhs.area - rhs.area, lhs.volume - rhs.volume, lhs.unit_bending_energy - rhs.unit_bending_energy);
  }

    //! Overloaded addition and assignment operator.
    /**
     * Works through the use of an already overloaded addition operator and
     * simply performs `lhs = lhs + rhs;`.
     *
     * @param lhs
     * @param rhs
     * @see Geometry<Real, Index> operator+(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
     */
    friend void operator+=(Geometry<Real, Index>& lhs, Geometry<Real, Index> const& rhs)
    {
        lhs = lhs + rhs;
    }

    //! Overloaded subtraction and assignment operator.
    /**
     * Works through the use of an already overloaded subtraction operator and
     * simply performs `lhs = lhs - rhs;`.
     *
     * @param lhs
     * @param rhs
     * @see Geometry<Real, Index> operator-(Geometry<Real, Index> const& lhs, Geometry<Real, Index> const& rhs)
     */
    friend void operator-=(Geometry<Real, Index>& lhs, Geometry<Real, Index> const& rhs){
        lhs = lhs - rhs;
    }

    //! Overloaded addition and assignment operator.
    /**
     * Can accumulate a node onto a Geometry struct.
     * @param node Geometric data from this node is added to the corresponding data members of the Geometry struct.
     */
    void operator+=(Node<Real, Index> const& node){
      area += node.area;
      volume += node.volume;
      unit_bending_energy += node.unit_bending_energy;
    }

};

/**
 * @GlobalsStub
 * @{
 */

//! This enum defines named types of triangulations that are implemented in flippy.
/**
 * A triangulation type needs to be provided as a template parameter to the Triangulation class during instantiation.
 * This will tell the class to create an appropriate topology of the triangulated Nodes.
 * The named types tell the triangulation class to...
 */
enum TriangulationType{

    //! Create a spherical triangulation which is a sub-triangulation of a regular icosahedron
    SPHERICAL_TRIANGULATION,
    //! Create a triangulation which is a sub-triangulation of a plane square.
    EXPERIMENTAL_PLANAR_TRIANGULATION
};
/**@}*/


/**
 * \brief Implementation of Triangulation of two-dimensional surfaces in 3D.
 * \image html assets/triangulation.png "Figure tr 1. Visualization of the triangulation."  width=500px
 * <b>Visualization of the triangulation</b>. <b>A</b>: Triangulated sphere with \f$N_{\mathrm{nodes}}=2252\f$. Black edges highlight the local neighborhood of a node.
 * Circular arrows show the counterclockwise orientation of the nodes.
 * This choice guarantees that all normal vectors point to the outside of the sphere.
 * B: An arbitrary node \f$i\f$, with its curvature vector \f$\vec{K}_{i}\f$ and a highlighted angle \f$\alpha^j_{i,j+1}\f$
 * at neighbour \f$j\f$ opposite to the edge \f$i, j+1\f$. Superscript \f$j\f$ denotes the neighboring node to which the angle
 * belongs and subscript \f$i,j+1\f$ denotes the edge opposite of the angle.
 * <b>C</b>: Node \f$i\f$ with its associated Voronoi area \f$A_i\f$ highlighted in red.
 * The node has an associated area inside each triangle it is part of.
 * We also highlight the triangle \f$i,j,j+1\f$ (light red with stripes) with the
 * face normal \f$\vec{n}_{i,j,j+1}\f$ and the area \f$A_{i,j,j+1}\f$.
 * The part of this triangle that is associated with node \f$i\f$ is highlighted in dark red and has the area \f$A_{ij}\f$.
 * The convention is to use the central and rightmost nodes in the subscript.
 * Since the nodes are ordered counterclockwise, this convention is unambiguous.
 * <b>D</b>: Volume associated with node \f$i\f$ is made up of tetrahedrons
 * that have as their base the triangles that make up the Voronoi cell of the node.
 * The head of the tetrahedron points to some lab frame origin \f$\cal{O}\f$. \f$V_{ij}\f$ is the part of the volume
 * associated to node \f$i\f$ that has its base in the triangle \f$i,j,j+1\f$.
 *
 *
 * @tparam Real @RealStub
 * @tparam Index @IndexStub
 * @tparam triangulation_type Template parameter that specifies the type of triangulation to be created.
 * This parameter must be chosen from the fp::TriangulationType `enum`.
 * Defaulted to fp::SPHERICAL_TRIANGULATION.
 */
template<floating_point_number Real, indexing_number Index, TriangulationType triangulation_type=SPHERICAL_TRIANGULATION>
class Triangulation
{
private:
    explicit Triangulation(Real verlet_radius_inp)
    :global_geometry_(), verlet_radius(verlet_radius_inp){}
public:
    Triangulation() = default;
    //unit tested
    //! Constructor that can re-initiate a triangulation from the stored data.
    /**
     *
     * @param nodes_input json object that contains data generated by make_egg_data() function, or similarly structured data.
     * @param verlet_radius_inp Value for [Verlet radius](https://en.wikipedia.org/wiki/Verlet_list).
     */
//    [[deprecated("This constructor is deprecated and will be removed in the future. Use Triangulation::load_sphere_from_json() instead.")]]
    Triangulation(Json const& nodes_input, Real verlet_radius_inp):Triangulation(verlet_radius_inp)
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

    //! Constructor that can initiate a spherical triangulation from scratch.
    /**
     *
     * @param n_nodes_iter Number of sub-triangulations.
     * @param R_initial_input Initial radius of the spherical triangulation.
     * @param verlet_radius_inp Value for [Verlet radius](https://en.wikipedia.org/wiki/Verlet_list).
     */
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

    //! Constructor that can initiate a planar triangulation from scratch.
    /**
     * This overload initiates a planar triangulation with free floating (non-periodic boundaries).
     * Edge nodes only serve a topological function and they do not contribute curvature or area to the membrane.
     * @warning @linearTriangulationWarningStub
     * @param n_length Number of nodes along the length of the triangulation (at the boundary).
     * @param n_width Number of nodes along the width of the triangulation (at the boundary).
     * @param length Length of the planar membrane.
     * @param width Width of the planar membrane
     * @param verlet_radius_inp Value for [Verlet radius](https://en.wikipedia.org/wiki/Verlet_list).
     */
    Triangulation(Index n_length, Index n_width, Real length, Real width, Real verlet_radius_inp):Triangulation(verlet_radius_inp)
    {
        static_assert(triangulation_type == EXPERIMENTAL_PLANAR_TRIANGULATION, "This initialization is intended for planar triangulations");
        triangulate_planar_nodes(n_length, n_width, length, width);
        orient_plane();
        initiate_advanced_geometry();
    }


    //! Special constructor that can initialize a triangulation with spherical topology from a binary stl [stl](https://en.wikipedia.org/wiki/STL_(file_format)#Format) file.
    /**
     * This constructor is intended for initializing a triangulation with spherical topology from an stl file.
     * @warning: This constructor is experimental and is not fully tested. The constructor trusts the stl file that it is a properly formed triangulation with spherical topology and a well oriented surface.
     * currently there are no catches in place to gracefully handle malformed stl files. malformed stl files will simply result in runtime errors.
     *
     * @param stl_file_path an absolute or relative path to a binary stl file
     * @param verlet_radius_inp
     * @return
     */
    static Triangulation<Real, Index, SPHERICAL_TRIANGULATION> experimental_load_sphere_from_stl(std::filesystem::path const& stl_file_path, Real verlet_radius_inp){
            std::vector<implementation::stlTriangle<Real, Index>> triangles = implementation::stlSerializer<Real, Index>::read_STLSolid_into_triangle_vec(stl_file_path);

            std::vector<Node<Real, Index>> nodes;
            std::vector<implementation::stlNode<Real, Index>> unique_nodes;
            Triangulation<Real, Index> triangulation;
            nodes.resize(triangles.size()/2 + 2);
            for (Index i=0; implementation::stlTriangle<Real, Index> &triangle : triangles) {
                for (unsigned int tr_idx = 0; tr_idx < 3; ++tr_idx) {
                    implementation::stlNode<Real, Index>& node = triangle[tr_idx];
                    auto node_pos = std::find(unique_nodes.begin(), unique_nodes.end(), node);
                    if (node_pos == unique_nodes.end()) {
                        node.id = i;
                        triangle[tr_idx].id = i;
                        ++i;
                        unique_nodes.push_back(node);
                        Node<Real, Index> n{};
                        n.pos = node.pos;
                        n.id = node.id;
                        nodes[n.id] = n;
                    }else{
                        triangle[tr_idx].id = node_pos->id;
                    }
                }
            }

            for (auto &triangle : triangles) {
                for (Index tr_idx = 0; tr_idx < 3; ++tr_idx) {
                    Index id_0 = triangle[tr_idx].id;
                    auto [id_1, id_2] = triangle.other_node_ids(tr_idx);

                    if(!is_member(nodes[id_0].nn_ids, id_1)){nodes[id_0].nn_ids.push_back(id_1);}
                    if(!is_member(nodes[id_0].nn_ids, id_2)){nodes[id_0].nn_ids.push_back(id_2);}

                }
            }
            Triangulation<Real, Index, SPHERICAL_TRIANGULATION> tr;
            tr.verlet_radius = verlet_radius_inp;
            tr.R_initial = 1;
            std::cout << "finished stl parsing\n";
            tr.nodes_ = Nodes<Real, Index>(nodes);
            tr.all_nodes_are_bulk();
            std::cout << "node count: " << tr.nodes_.size() << "\n";
            tr.orient_surface_of_a_sphere();
            std::cout << "surface oriented\n";
            tr.initiate_advanced_geometry();
            std::cout << "advanced geometry initiated\n";
            return tr;
        }


    //! Set the radius of the Verlet list to a new value.
    /**
     * @param R new radius
     * @related Triangulation::make_verlet_list()
     */
    void set_verlet_radius(Real R){
        verlet_radius = R;
        verlet_radius_squared = R*R;
    }

    //todo unittest
    //! Create a [Verlet list](https://en.wikipedia.org/wiki/Verlet_list).
    /**
     * This method creates a Verlet list for each node of the triangulation. All nodes that are inside the `verlet_radius`, of a given node,
     * are included in its Verlet list (Node.verlet_list),
     */
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

    //! Adds the same 3D vector to the positions of each node of the triangulation.
    /**
     * This method is most helpful in shifting a triangulation after its initiation.
     * For example, to set up initial conditions for a simulation where several triangulated vesicles interact.
     * @param translation_vector A fixed 3D vector by which the triangulation is to be shifted.
     */
    void translate_all_nodes(vec3<Real> const& translation_vector)
    {
        for (Index i = 0; i<nodes_.size(); ++i) { move_node(i, translation_vector); }
    }

    //unit tested
    //! Calculate the area-weighted average of node positions.
    /**
     * This method calculates the average position of the mass center of the triangulation,
     * by averaging the positions of each node of the triangulation.
     * The average is weighted by the area associated with each node.
     * @return position of the mass center.
     */
    vec3<Real> calculate_mass_center() const
    {
        vec3<Real> mass_center = vec3<Real>{0., 0., 0.};
        for (auto const& node : nodes_) { mass_center += node.pos; }
        mass_center = mass_center/static_cast<Real>(nodes_.size());
        return mass_center;
    }

    //unit tested
    //! Move an individual node of the triangulation and update all the geometric quantities of the triangulation that changed.
    /**
     * @param node_id @NodeIDStub
     * @param displacement_vector 3D vector by which the chosen node is to be displaced.
     */
    void move_node(Index node_id, vec3<Real> const& displacement_vector)
    {
        pre_update_geometry = get_two_ring_geometry(node_id);
        nodes_.displace(node_id, displacement_vector);
        update_two_ring_geometry(node_id);
        post_update_geometry = get_two_ring_geometry(node_id);
        update_global_geometry(pre_update_geometry, post_update_geometry);
    }

    // unit-tested
    //! Adds a new node to the next neighbor list of a given node and calculates their mutual distance.
    /**
     * @note The direct use of this method is discouraged unless a very specific type of edge flipping is required.
     * Triangulation::flip_bond is a more high-level method, suitable for most basic bond-flipping needs.
     *
     * This method finds the anchor node in the Nodes::nn_ids vector of the center_node
     * and uses Node classes own method Node::emplace_nn_id to emplace the new_value
     * there (together with its distance to the center_node).
     * @param center_node_id @NodeIDStub The new node is emplaced in the Node.nn_ids vector of this node.
     * @param anchor_id @LocNNIndexStub This is the index of the next neighbor (inside `nn_ids` vector of the node `center_node_id`), before which the new node id is emplaced.
     * @param new_value @NodeIDStub The id of the new node.
     */
    void emplace_before(Index center_node_id, Index anchor_id, Index new_value)
    {
        // The body of this function looks like it does not guard against find returning
        // end() pointer, but this is taken care of in the emplace_nn_id method.
        auto anchor_pos_ptr = std::find(nodes_[center_node_id].nn_ids.begin(),
                nodes_[center_node_id].nn_ids.end(), anchor_id);
        indexing_number auto anchor_pos = (Index) (anchor_pos_ptr - nodes_[center_node_id].nn_ids.begin());
        nodes_[center_node_id].emplace_nn_id(new_value, nodes_[new_value].pos, anchor_pos);
    }

    /** \brief Securely flip the bond inside a quadrilateral formed by the nodes given by node_id,
     * nn_id and their two common next neighbors, if all topological requirements are satisfied.
     */
    BondFlipData<Index> flip_bond(Index node_id, Index nn_id,
                                  Real min_bond_length_square,
                                  Real max_bond_length_square){
   /**
     * *flip_bond* function takes a lot of care to keep the triangulation intact, i.e., not introduce holes or additional bonds in it.
     * To this end, the function performs a lot of checks to determine if the proposed bond fip is allowed by the topology.
     * The information if the flip was allowed and succeeded or not will be encoded in the BondFlipData struct, which this function will return.
     * The struct also contains information on the new end nodes of the bond, which is useful if one wishes to undo the flip.
     * A bond flip can fail if:
     * - if the provided nodes id's do not correspond to neighboring nodes.
     * - one of the donor nodes have already too few bonds (bonds less than #BOND_DONATION_CUTOFF).
     * - if a new bond that would be created, would be
     *  - too long (length squared is larger than **max_bond_length_square**)
     *  - too short (length squared is smaller than **min_bond_length_square**)
     *
     * @param node_id @NodeIDStub
     * @param nn_id @NNIDStub
     * @param min_bond_length_square @BondLengthSquareStub{minimal}
     * @param max_bond_length_square @BondLengthSquareStub{maximal}
     * @return If the flip was successful, the returned BondFlipData struct will contain a boolean state variable BondFlipData::flipped = **true**.
     * The members BondFlipData::common_nn_0 and BondFlipData::common_nn_1 will contain global ids of new end nodes of the flipped bond.
     * If the flip was not successful, then a default initialized BondFlipData struct will be returned with BondFlipData::flipped = **false**.
     * @note Regardless of the return values, the primary purpose of the function, that of flipping a bond, is accomplished as a side-effect.
     */
        if constexpr (triangulation_type == TriangulationType::SPHERICAL_TRIANGULATION) {
            return flip_bulk_bond(node_id, nn_id, min_bond_length_square, max_bond_length_square);
        } else if constexpr (triangulation_type == TriangulationType::EXPERIMENTAL_PLANAR_TRIANGULATION){

            if (boundary_nodes_ids_set_.contains(node_id) or boundary_nodes_ids_set_.contains(nn_id)){
            }else{
                Neighbors<Index> common_nns = previous_and_next_neighbour_global_ids(node_id, nn_id);
                if(boundary_nodes_ids_set_.contains(common_nns.j_m_1) or boundary_nodes_ids_set_.contains(common_nns.j_p_1)){
                }else{
                    return flip_bond_in_quadrilateral(node_id, nn_id, common_nns, min_bond_length_square, max_bond_length_square);
                }
            }
            return BondFlipData<Index>();

        }else{
            static_assert(triangulation_type == TriangulationType::SPHERICAL_TRIANGULATION or triangulation_type == TriangulationType::EXPERIMENTAL_PLANAR_TRIANGULATION,
                          "Triangulation type must be either spherical or planar.");
        }
    }

    //unit-tested
    //! Un-flip a bond that was just flipped.
    /**
     * This method reverses a flip between two nodes that used to be connected and had their node just flipped away,
     * provided their id's and the BondFlipData that holds the information on the current bond holders.
     *
     *
     * This means that `node_id` and `nn_id` are (pre-flip) owners of the bond. And `common_nns` contains the id's of the
     * current bond owners (to which the bond was flipped to).
     * The edge will be taken away from the current two neighbors and added to the previous owners. All geometric information of the triangulation will be updated.
     *
     * @warning This method only works correctly if there were no bond flips after that flip which is being reversed.
     * I.e., this method can only reverse the last flip. Furthermore, the method assumes that the arguments are provided correctly and
     * will result in an illegal triangulation if this is not the case.
     *
     * @note This method does not check the validity of its input to provide a fast way of flip reversal. It relies on the fact that
     * the user usually has the exact knowledge of all four node id's that participated in a flip and can provide those id's
     * in the correct order. For a safer way of un-fliping, the Triangulation::flip_bond method can be used twice, with interchanged order
     * of arguments.
     *
     * The following example implementation of a Monte Carlo flip update method shows the proper way to use
     * the Triangulation.flip_bond and Triangulation.unflip_bond methods.
     *
     * ```c++
     *  // `e_old`, `e_new`, `triangulation and `parameters` are data members of the updater
     *  void mc_flipp_update(fp::Node<Real, Index> const& node)
     *  {
     *       e_old = energy_function(node, triangulation, parameters);
     *       Index number_nn_ids = node.nn_ids.size();
     *       Index nn_id = node.nn_ids[std::uniform_int_distribution<Index>(0, number_nn_ids-1)(rng)];
     *       auto bond_flip_data = triangulation.flip_bond(node.id, nn_id, min_bond_length_square, max_bond_length_square);
     *       if (bond_flip_data.flipped) {
     *           e_new = energy_function(node, triangulation, parameters);
     *           if (move_needs_undoing()) { triangulation.unflip_bond(node.id, nn_id, bond_flip_data); }
     *       }
     *   }
     *```
     *@see Triangulation.flip_bond BondFlipData
     *
     * @param node_id @NodeIDStub
     * @param nn_id @NNIDStub
     * @param common_nns BondFlipData containing information on the common next neighbor ids.
     */
    void unflip_bond(Index node_id, Index nn_id, BondFlipData<Index> const& common_nns)
    {
        flip_bond_unchecked(common_nns.common_nn_0, common_nns.common_nn_1, nn_id, node_id);
        update_diamond_geometry(node_id, nn_id, common_nns.common_nn_0, common_nns.common_nn_1);
        update_global_geometry(post_update_geometry, pre_update_geometry);
    }

    //! Exchange the next neighborhood between four nodes in a manner that will correspond to
    //! a bond flip if the provided information was correct.
    /**
     * @warning This method may lead to an unphysical state (broken lattice that no longer represents a triangulation) if the provided
     * arguments are not correct.
     * @note For most use-cases the high level methods Triangulation::flip_bond and Triangulation::unflilp_bond are recommended.
     * They guarantee that the performed flips result in a legal triangulation and reject the flip otherwise.
     * This method should only be used by advanced users that have very specific needs for a flip updater.
     *
     * The correct functioning of this method requires that, `node_id` and `nn_id` need to be next neighbours, and `common_nn_j_m_1` and `common_nn_j_p_1`
     * need to be common next neighbors of both `node_id` and `nn_id`.
     * Moreover, the nodes need to be ordered in the Node::nn_ids vector of the Node represented by `node_id` as follows:
     * `... common_nn_j_m_1, node_id, common_nn_j_p_1 ...`,
     * or a cyclic permutation thereof.
     * @param node_id @NodeIDStub
     * @param nn_id @NNIDStub
     * @param common_nn_j_m_1
     * @param common_nn_j_p_1
     * @return BondFlipData, where the flipped field is always set to `True`.
     */
    BondFlipData<Index> flip_bond_unchecked(Index node_id, Index nn_id,
                                            Index common_nn_j_m_1, Index common_nn_j_p_1)
    {
        emplace_before(common_nn_j_m_1, node_id, common_nn_j_p_1);
        emplace_before(common_nn_j_p_1, nn_id, common_nn_j_m_1);
        delete_connection_between_nodes_of_old_edge(node_id, nn_id);
        return {.flipped=true, .common_nn_0=common_nn_j_m_1, .common_nn_1=common_nn_j_p_1};
    }

    // unit-tested
    //! Update the geometric quantities associated with the given node.
    /**
     * This is the core update function of the triangulation.
     * It updates the geometric quantities of the given node after its position in the triangulation has been changed
     * through a node move, or its neighborhood structure has been changed through a bond flip.
     * Calculate and update the local curvature, area, volume, and unit bending energy of the node (See Figure tr1 B, C and D).
     *
     * @see Node::curvature_vec, Node::area, Node::volume, Node::unit_bending_energy
     * @param node_id @NodeIDStub
     */
    void update_bulk_node_geometry(Index node_id)
    {
        update_nn_distance_vectors(node_id);

        Real area_sum = 0.;
        vec3<Real> face_normal_sum{0., 0., 0.}, local_curvature_vec{0., 0., 0.};
        vec3<Real> face_normal;
        indexing_number auto nn_number = (Index) nodes_.nn_ids(node_id).size();
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
#ifdef DEBUG
            if(face_normal_norm < 1e-10) {
                throw std::runtime_error("A triangle face is degenerate and Area sum is evaluating to "+std::to_string(face_normal_norm)+". This should not happen.");
            }
#endif
            face_area = mixed_area(lij, lij_p_1, Real(0.5)*face_normal_norm, cot_at_j, cot_at_j_p_1);
            area_sum += face_area;
            face_normal_sum += face_area*face_normal/face_normal_norm;

            local_curvature_vec -= (cot_at_j_p_1*lij + cot_at_j*lij_p_1);
        }
        nodes_.set_area(node_id, area_sum);
        nodes_.set_volume(node_id, nodes_[node_id].pos.dot(face_normal_sum)/((Real) 3.)); // 18=3*6: 6 has the aforementioned justification. 3 is part of the formula for the tetrahedron volume
        nodes_.set_curvature_vec(node_id,  -local_curvature_vec/((Real) 2.*area_sum)); // 2 is part of the formula to calculate the local curvature I just did not divide the vector inside the loop
        nodes_.set_unit_bending_energy(node_id, local_curvature_vec.dot(local_curvature_vec)/((Real) 8.*area_sum)); // 8 is 2*4, where 4 is the square of the above two and the area in the denominator is what remains after canceling. 1/ comes from the pre-factor to bending energy

    };


    //! This function is deprecated!
    /**
     * @warning This function is here for legacy reasons! It is deprecated and will be removed in future minor updates.
     * This function returns values of eqn.'s (82) & (84) from the paper [Gueguen et al. 2017](https://doi.org/10.1039/C7SM01272A).
     * Every node has its associated Voronoi area, and each Voronoi area can be subdivided into parts that are
     * associated with each triangle that the node is part of. This function returns that sub-area and the face normal
     * of that triangle.
     * @param lij distance between a node and its jth next neighbor
     * @param lij_p_1 distance between a node and its j+1th next neighbor
     * @return tuple of the area associated with the node inside the triangle (i,j,j+1) and the face normal of the triangle.
     */
    [[deprecated("This function is deprecated and will be removed in a future release. It uses expensive function calls and is not recommended for use.")]]
    static std::tuple<Real, vec3<Real>> partial_voronoi_area_and_face_normal_of_node_in_a_triangle(vec3<Real> const& lij,
                                                                                                   vec3<Real> const& lij_p_1)
    {
        Real area, face_normal_norm;
        vec3<Real> un_noremd_face_normal;
        //precalculating this normal and its norm will be needed in area calc. If all triangles are oriented as
        // right-handed, then this normal will point outwards
        un_noremd_face_normal = lij.cross(lij_p_1);
        face_normal_norm = un_noremd_face_normal.norm();
        area = mixed_area(lij, lij_p_1, face_normal_norm/2.);
        return std::make_tuple(area, un_noremd_face_normal);
    }

    //! The node-associated area inside a triangle.
    /**
     * This function is calculating the area associated with a node inside a triangle. As depicted in Figure tr1. A.
     * This can be found in the description of the Triangulation class.
     *
     * Every node of the triangulation has its own associated area.
     * In the simplest case, the area associated with a node is the area of the Voronoi cell of that node. Where the Voronoi tessellation of the surface is the dual lattice of the triangulation.
     * This function returns the area associated with a node inside a triangle. If that triangle is not obtuse, then this area is simply
     * the area of the Voronoi cell of the node that is inside the triangle.
     * This Voronoi area becomes negative for obtuse triangles, and thus an exception has to be made.
     * If the triangle has an obtuse angle at the node, then the area associated with the node inside the triangle is half of the triangle area,
     * otherwise, it is the quarter of the triangle area. This procedure is described in detail by [Meyer et al. 2003](https://doi.org/10.1007/978-3-662-05105-4_2)
     * @param lij Distance vector between the node and its next neighbor. (Next neighbors are ordered according to the right-hand rule. See Figure tr1. A and C)
     * @param lij_p_1 Distance vector between the node and its next neighbor. (Next neighbors are ordered according to the right-hand rule)
     * @param triangle_area area of the triangle i,j,j+1. See Area \f$A_{i,j,j+1}\f$ in Figure tr1. C.
     * @param cot_at_j Cotangent of the angle at the node j, opposite to the edge i,j+1. See Figure tr1. B.
     * @param cot_at_j_p_1 Cotangent of the angle at the node j+1, opposite to the edge i,j. See Figure tr1. B.
     * @return Area associated with the node inside the triangle. See Area \f$A_{ij}\f$ in Figure tr1. C.
     */
    static Real mixed_area(vec3<Real> const& lij, vec3<Real> const& lij_p_1, Real triangle_area, Real cot_at_j, Real cot_at_j_p_1){
        if ((cot_at_j>0.) && (cot_at_j_p_1>0.)) { // both angles at j and j+1 are smaller than 90 deg so the triangle can only be obtuse at the node
            if (lij.dot(lij_p_1)>0) { // cos at i is positive i.e. angle at i is not obtuse
                return (cot_at_j_p_1*lij.dot(lij) + cot_at_j*lij_p_1.dot(lij_p_1))/Real(8.);
            }
            else {//obtuse at node i.
                return triangle_area/Real(2.);
            }
        }
        else {//obtuse at node j or j+1.
            return triangle_area/Real(4.);
        }

        }

    //unit tested
    //! This function is deprecated!
    /**
     * @warning This function is here for legacy reasons! It is deprecated and will be removed in future minor updates.
     * Use the alternative mixed_area function instead!
     *
     * @param lij Distance vector between the node and its next neighbor. (Next neighbors are ordered according to the right-hand rule. See Figure tr1. A and C)
     * @param lij_p_1 Distance vector between the node and its next neighbor. (Next neighbors are ordered according to the right-hand rule)
     * @param triangle_area area of the triangle \f$i,j,j+1\f$. See Area \f$A_{i,j,j+1}\f$ in Figure tr1. C.
     * @return Area associated with the node inside the triangle. See Area \f$A_{ij}\f$ in Figure tr1. C.
     */
    [[deprecated("This function is deprecated and will be removed in a future release. mixed_area, which does not take precalculated cotangents, performs expensive calculations! Use the alternative mixed_area function!")]]
    static Real mixed_area(vec3<Real> const& lij, vec3<Real> const& lij_p_1, Real const& triangle_area)
    {
        vec3<Real> ljj_p_1 = lij_p_1 - lij;

        Real cot_at_j = cot_between_vectors(lij, (-1)*ljj_p_1);
        Real cot_at_j_p_1 = cot_between_vectors(lij_p_1, ljj_p_1);
        if ((cot_at_j>Real(0.)) && (cot_at_j_p_1>Real(0.))) { // both angles at j and j+1 are smaller than 90 deg so the triangle can only be obtuse at the node
            if (lij.dot(lij_p_1)>Real(0.)) { // cos at i is positive i.e. angle at i is not obtuse
                return (cot_at_j_p_1*lij.dot(lij) + cot_at_j*lij_p_1.dot(lij_p_1))/8.;
            }
            else {//obtuse at node i.
                return triangle_area/Real(2.);
            }
        }
        else {//obtuse at node j or j+1.
            return triangle_area/Real(4.);
        }

    }

    //! Aggregates and Returns the geometric quantities of the center node and its next neighbor nodes.
    /**
     * This function aggregates area volume and squared curvature integrated over the area for the  two-ring of the
     * associated with the node, using the stored quantities.
     * Two-ring refers to the second concentric ring of the next-nearest-neighbor nodes surrounding the center node.
     *
     * @param node_id
     * @return Geometry<Real, Index> object containing the geometric quantities of the center node and its next neighbor nodes.
     */
    [[nodiscard]] Geometry<Real, Index> get_two_ring_geometry(Index node_id) const
    {
        Geometry<Real, Index> trg(nodes_[node_id]);
        for (auto const& nn_id: nodes_[node_id].nn_ids) {
            trg += nodes_[nn_id];
        }
        return trg;
    }

    //! Updates the geometric quantities of the center node and its next neighbor nodes.
    /**
     * This function calculates and updates the area volume and squared curvature integrated over the area for the  two-ring of the
     * associated with the node.
     * Two-ring refers to the second concentric ring of the next-nearest-neighbor nodes surrounding the center node.
     * @param node_id @NodeIDStub
     */
    void update_two_ring_geometry(Index node_id)
    {
        if constexpr(triangulation_type == TriangulationType::SPHERICAL_TRIANGULATION) {
            update_two_ring_geometry_on_a_boundary_free_triangulation(node_id);
        }
        else if constexpr(triangulation_type == TriangulationType::EXPERIMENTAL_PLANAR_TRIANGULATION) {
            update_two_ring_geometry_on_a_boundary_triangulation(node_id);
        }
        else {
            static_assert(triangulation_type == TriangulationType::SPHERICAL_TRIANGULATION || triangulation_type == TriangulationType::EXPERIMENTAL_PLANAR_TRIANGULATION,
                          "Triangulation type is not supported!");
        }
    };

    // unit-tested
    //! Method for stretching or squeezing the initial triangulation shape.
    /**
     * This method is most useful for transforming a spherical triangulation into an ellipse.
     *
     * The method stretches `x`, `y`, and `z` components of each node by a factor provided in the function argument.
     * Triangulation's local and global geometric properties are updated after the stretch.
     *
     * @note The stretch will happen with respect to the lab frame and not the mass center of the triangulation.
     * Thus, to get the intended results, the stretch is most likely desired,
     * when the triangulation is centered around the origin of the lab frame, This is the case right after the initiation.
     *
     * @param x_stretch Stretching factor of the `x` component of the position of the triangulation nodes.
     * @param y_stretch Stretching factor of the `y` component of the position of the triangulation nodes.
     * @param z_stretch Stretching factor of the `z` component of the position of the triangulation nodes.
     */
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
    //! Aggregates the geometric quantities of the diamond configuration of nodes associated with a bond flip.
    /**
     *```txt
     * State of the diamond, which contains all nodes that participate in a bond flip.
     *  before the flip
     *
     *      common nn 1
     *      /          \
     *     /            \
     *   node --------- nn
     *     \            /
     *      \          /
     *     common nn 0
     *```
     *
     * @param node_id @NodeIDStub
     * @param nn_id @NodeIDStub
     * @param cnn_0 Global id of the common next nearest neighbor of node_id and nn_id.
     * @param cnn_1 Global id of the common next nearest neighbor of node_id and nn_id.
     * @return Geometric quantities aggregated over the diamond configuration of nodes associated with a bond flip.
     */
    [[nodiscard]] Geometry<Real, Index> calculate_diamond_geometry(Index node_id, Index nn_id,
                                                             Index cnn_0, Index cnn_1) const
    {
        Geometry<Real, Index> diamond_geometry(nodes_[node_id]);
        diamond_geometry += nodes_[nn_id];
        diamond_geometry += nodes_[cnn_0];
        diamond_geometry += nodes_[cnn_1];
        return diamond_geometry;
    };

    //Todo unittest
    //! Calculates and updates the geometric quantities of the diamond configuration of nodes associated with a bond flip.
    /**
     *```txt
     * State of the diamond, which contains all nodes that participate in a bond flip.
     *  before the flip
     *
     *      common nn 1
     *      /          \
     *     /            \
     *   node --------- nn
     *     \            /
     *      \          /
     *     common nn 0
     *```
     *
     * @param node_id @NodeIDStub
     * @param nn_id @NodeIDStub
     * @param cnn_0 Global id of the common next nearest neighbor of node_id and nn_id.
     * @param cnn_1 Global id of the common next nearest neighbor of node_id and nn_id.
     */
    void update_diamond_geometry(Index node_id, Index nn_id, Index cnn_0, Index cnn_1)
    {
        update_bulk_node_geometry(node_id);
        update_bulk_node_geometry(nn_id);
        update_bulk_node_geometry(cnn_0);
        update_bulk_node_geometry(cnn_1);
    };

    // Const Viewer Functions
    //! Returns the number of nodes in the triangulation.
    /**
     * @return Number of the nodes in the triangulation.
     */
    [[nodiscard]] Index size() const { return nodes_.size(); }
    //! Returns a constant reference to the node with the given id.
    /**
     * Triangulation will never give non-constant access to a node.
     * In order to change a node, one has to use the methods of the Triangulation class.
     * This guarantees that the triangulation is always in a consistent state.
     * @param idx @NodeIDStub
     * @return Constant reference to the node with the given id.
     */
    const Node<Real, Index>& operator[](Index idx) const { return nodes_.data.at(idx); }
    //! Returns a constant reference to the underlying Nodes container.
    /**
     * @return Constant reference to the underlying Nodes container.
     */
    const Nodes<Real, Index>& nodes() const { return nodes_; }
    //! Creates a JSON object with the data of the triangulation.
    /**
     * Egg refers to the fact that the data can be used to recreate the triangulation using the Triangulation(Json const& nodes_input, Real verlet_radius_inp) constructor.
     * @note The Triangulation(Json const& nodes_input, Real verlet_radius_inp) constructor is currently only implemented for a spherical Triangulation!
     *
     * @return Triangulation data in JSON format.
     */
    [[nodiscard]] Json make_egg_data() const { return nodes_.make_data(); }
    //! Information about the global geometric quantities of the triangulation, like global area, volume, and total unit bending energy.
    /**
     * @return Geometric quantities of the triangulation aggregated over all nodes.
     */
    [[nodiscard]] const Geometry<Real, Index>& global_geometry() const { return global_geometry_; }

    //Todo unittest
    //! Initiates the global geometry of the triangulation.
    /**
     * The global geometry is calculated by summing up the local geometries of all nodes.
     *
     */
    void make_global_geometry()
    {
        const Geometry<Real, Index> empty{};
        global_geometry_ = empty;
        for (auto node_id: bulk_nodes_ids) {
            update_bulk_node_geometry(node_id);
            update_global_geometry(empty, Geometry<Real, Index>(nodes_[node_id]));
        }
        for (auto node_id: boundary_nodes_ids_set_) {
            update_boundary_node_geometry(node_id);
            update_global_geometry(empty, Geometry<Real, Index>(nodes_[node_id]));
        }
    }

    //! Returns the ids of all nodes that are not on the boundary.
    /**
     * Only works for triangulations that have a boundary.
     * @return unique set of global ids of all nodes that are not on the boundary.
     */
    std::set<Index> boundary_nodes_ids_set() const {
        static_assert(triangulation_type == TriangulationType::EXPERIMENTAL_PLANAR_TRIANGULATION, "This function is only implemented for PLANAR_TRIANGULATION.");
        return boundary_nodes_ids_set_;
    }

    //Todo unittest
    //! Updates the local geometry of a boundary node (for triangulation types that have a boundary).
    /**
     * @warning @linearTriangulationWarningStub
     * Boundary nodes need to be treated differently, depending on the boundary conditions.
     * Right now, flippy is only handling fixed boundary conditions, where all geometric quantities of boundary nodes are set to zero.
     * @param node_id Id of a boundary node.
     */
    void update_boundary_node_geometry([[maybe_unused]]Index node_id){
        update_nn_distance_vectors(node_id);
    }

#ifdef TESTING_FLIPPY_TRIANGULATION_ndh6jclc0qnp274b
public:
#else
private:
#endif
    Real R_initial;
    Nodes<Real, Index> nodes_;
    std::vector<Index> bulk_nodes_ids;
    Geometry<Real, Index> global_geometry_;
    Geometry<Real, Index> pre_update_geometry, post_update_geometry;
    mutable vec3<Real> l0_, l1_;
    Real verlet_radius{};
    Real verlet_radius_squared{};
    std::set<Index> boundary_nodes_ids_set_;

    //unit tested
    void initiate_advanced_geometry(){
        initiate_distance_vectors();
        make_global_geometry();
        set_verlet_radius(verlet_radius);
        make_verlet_list();
    }

    void update_two_ring_geometry_on_a_boundary_free_triangulation(Index node_id){
        update_bulk_node_geometry(node_id);
        for (auto nn_id: nodes_.nn_ids(node_id)) {
            update_bulk_node_geometry(nn_id);
        }
    }
    void update_two_ring_geometry_on_a_boundary_triangulation(Index node_id){
        if(boundary_nodes_ids_set_.contains(node_id)){
            update_boundary_node_geometry(node_id);
        }else{
            update_bulk_node_geometry(node_id);
        }

        for (auto nn_id: nodes_.nn_ids(node_id)) {
            if(boundary_nodes_ids_set_.contains(nn_id)){
                update_boundary_node_geometry(nn_id);
            }else{
                update_bulk_node_geometry(nn_id);
            }
        }
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

    //! This function calculates distance vectors from a node to all of its neighbors.
    void update_nn_distance_vectors(Index node_id)
    {
        /**
         *  The directions of the distance vectors are (radiating outward) pointing from the node to neighbors.
         *  The function also preserves the order of neighbors. Meaning that the order of distance vectors is in the same
         *  order as the provided list of neighbor ids.
         *  @param node_id Global id of a node.
         */

        for (Index i = 0; auto nn_id: nodes_.nn_ids(node_id)) {
            nodes_.set_nn_distance(node_id, i, nodes_.pos(nn_id) - nodes_.pos(node_id));
            ++i;
        }
    }

    /**
     * given a node `i` and its neighbor `j`, they will share two common neighbor nodes, `p` and `m`.
     * This function finds the angles at `p` & `m` opposite of `i-j` link.
     * This function implements the cot(alpha_ij) + cot(beta_ij) from fig. (6c) from [.
     * The order of these neighbors does not matter for the correct sign of the angles.
     * @param node_id @NodeIDStub
     * @param nn_id @NNIDStub
     * @param cnn_0 common neighbor node 0
     * @param cnn_1 common neighbor node 1
     * @return cot(alpha_ij_jm1) + cot(alpha_ij_jp1)
     */
    Real cot_alphas_sum(Index node_id, Index nn_id, Index cnn_0, Index cnn_1) const
    {
        /**
         *
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

    void orient_plane()
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
            static_assert(triangulation_type == EXPERIMENTAL_PLANAR_TRIANGULATION, "This function is only well defined for a planar triangulation");
            std::vector<Index> nn_ids_temp;
            vec3<Real> li0, li1;
            vec3<Real> mass_center = calculate_mass_center();
            mass_center.z+=10;
            for (Index node_id = 0; node_id < nodes_.size(); ++node_id) { //ToDo modernize this loop
                if(boundary_nodes_ids_set_.contains(node_id)){
                    continue;
                }
                nn_ids_temp = order_nn_ids(node_id);
                li0 = nodes_[nn_ids_temp[0]].pos - nodes_[node_id].pos;
                li1 = nodes_[nn_ids_temp[1]].pos - nodes_[node_id].pos;
                if ((li0.cross(li1)).dot(nodes_[node_id].pos - mass_center) < 0) {
                    std::reverse(nn_ids_temp.begin(), nn_ids_temp.end());
                }
                nodes_.set_nn_ids(node_id, nn_ids_temp);
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
        static const Index vln = static_cast<const Index>(VERY_LARGE_NUMBER_);
        std::array<Index, 2> res{vln, vln};
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
        indexing_number auto nn_number = (Index)nodes_.nn_ids(node_id_0).size();
        Index j_p_1 = Neighbors<Index>::plus_one(j, nn_number);
        Index j_m_1 = Neighbors<Index>::plus_one(j, nn_number);
        std::array<Index, 2> res{nodes_.nn_id(node_id_0,j_m_1),
                nodes_.nn_id(node_id_0,j_p_1)};
        return res;
    }

    std::array<Index, 2> two_common_neighbour_positions(Index node_id_0, Index node_id_1) const
    {
        static const Index vln = static_cast<const Index>(VERY_LARGE_NUMBER_);
        std::array<Index, 2> res{vln, vln};
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
         *     	given i and j, this function finds the local ids of j-1 and j+1 nodes and returns them IN THAT ORDER;
         *     	This function relies on the fact that i & j are neighbors and will throw a nasty runtime error if they are
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
         *     	given i and j, this function finds the global ids of j-1 and j+1 nodes and returns them IN THAT ORDER;
         *     	This function relies on the fact that i & j are neighbors and will throw a nasty runtime error if they are
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

        fp::indexing_number auto nNewNodesOnEdge = static_cast<Index>(n_iter - 1);
        fp::indexing_number auto nBulk = static_cast<Index>(nNewNodesOnEdge*(nNewNodesOnEdge+1)/2);
        fp::indexing_number auto nNodes = static_cast<Index>(fp::implementation::IcosahedronSubTriangulation<Real,Index>::N_ICOSA_NODEs
            + fp::implementation::IcosahedronSubTriangulation<Real,Index>::N_ICOSA_EDGEs*n_iter
            + fp::implementation::IcosahedronSubTriangulation<Real,Index>::N_ICOSA_FACEs*nBulk);
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
        Node<Real, Index> node{};
        node.curvature_vec=fp::vec3<Real>{0.,0.,0.};
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
            else{boundary_nodes_ids_set_.insert(node_id);}
        }
    }
    void all_nodes_are_bulk(){
        for(auto const& node: nodes_){
            bulk_nodes_ids.push_back(node.id);
        }
    }


        //unit tested
        BondFlipData<Index> flip_bulk_bond(Index node_id, Index nn_id,
                                           Real min_bond_length_square,
                                           Real max_bond_length_square) {
            BondFlipData<Index> bfd{};
            if (nodes_.nn_ids(node_id).size() > BOND_DONATION_CUTOFF) {
                if (nodes_.nn_ids(nn_id).size() > BOND_DONATION_CUTOFF) {
                    Neighbors<Index> common_nns = previous_and_next_neighbour_global_ids(node_id, nn_id);
                    Real bond_length_square = (nodes_.pos(common_nns.j_m_1) - nodes_.pos(common_nns.j_p_1)).norm_square();
                    if ((bond_length_square < max_bond_length_square) && (bond_length_square > min_bond_length_square)) {
                        if (common_neighbours(node_id, nn_id).size() == 2) {
                            pre_update_geometry = calculate_diamond_geometry(node_id, nn_id, common_nns.j_m_1,
                                                                             common_nns.j_p_1);
                            bfd = flip_bond_unchecked(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                            if (common_neighbours(bfd.common_nn_0, bfd.common_nn_1).size() == 2) {
                                update_diamond_geometry(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                                post_update_geometry = calculate_diamond_geometry(node_id, nn_id, common_nns.j_m_1,
                                                                                  common_nns.j_p_1);
                                update_global_geometry(pre_update_geometry, post_update_geometry);
                            } else {
                                flip_bond_unchecked(bfd.common_nn_0, bfd.common_nn_1, nn_id, node_id);
                                bfd.flipped = false;
                            }
                        }
                    }
                }
            }
            return bfd;
        }

    ///// EXPERIMENTAL SUPPORT /////
    BondFlipData <Index>
    flip_bond_in_quadrilateral(Index node_id, Index nn_id, const Neighbors <Index> &common_nns,
                               Real min_bond_length_square, Real max_bond_length_square) {
        BondFlipData<Index> bfd{};
        Real bond_length_square = (nodes_.pos(common_nns.j_m_1) - nodes_.pos(common_nns.j_p_1)).norm_square();
        if ((bond_length_square<max_bond_length_square) && (bond_length_square>min_bond_length_square)) {
            if (common_neighbours(node_id, nn_id).size() == 2) {
                pre_update_geometry = calculate_diamond_geometry(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                bfd = flip_bond_unchecked(node_id, nn_id, common_nns.j_m_1, common_nns.j_p_1);
                if (common_neighbours(bfd.common_nn_0, bfd.common_nn_1).size() == 2) {
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
        return bfd;
    }

};

}
#endif //FLIPPY_TRIANGULATION_HPP

