#ifndef FLIPPY_NODES_EXPERIMENTAL_HPP
#define FLIPPY_NODES_EXPERIMENTAL_HPP
/**
 * @file
 * @brief This file contains the fp::Node and fp::Nodes classes, data structures that represent a single node of the
 * triangulation and the collection of all nodes of the triangulation, respectively.
 */
#include "external/json.hpp"
#include "vec3.hpp"
#include <unordered_set>
#include <vector>

namespace fp::experimental {

constexpr unsigned int MAX_ALLOWED_NEIGHBORS = 16;

template <typename T> class NeighborsArray : public std::array<T, MAX_ALLOWED_NEIGHBORS> {
    public:
    using std::array<T, MAX_ALLOWED_NEIGHBORS>::array;
    size_t size() { return size_; }

    void emplace(size_t idx, T &&val) {
        for (size_t i = size_; i > idx; --i) {
            this->data()[i] = this->data()[i - 1];
        }
        this->data()[idx] = val;
        size_++;
    }
    //    void pop(size_t idx){
    //        for(size_t i = size_ - 1; )
    //        move_left(idx, this);
    //    }

    private:
    //    void move_left(size_t start, )
    size_t size_{};
};

using Json = nlohmann::json;
//! A data structure containing all geometric and topological information associated with a node.
/**
 * This is a DUMB DATA STRUCTURE, meaning that it is not responsible for the coherence of the data it contains.
 * For performance reasons, methods associated with Node struct will never check if the Node::curvature is the norm of
 * the Node::curvature_vector or if the Node::nn_ids and Node::nn_distances are in the correct order. It is the
 * responsibility of higher-order structures like Nodes and Triangulation to check that correct data is stored and
 * updated correctly. However, it does check the data for consistency. It will match the length of Node::nn_ids and
 * Node::nn_distances and pop and add both of them together.
 * @tparam Real @RealStub
 * @tparam Index @IndexStub
 */
template <floating_point_number Real, indexing_number Index> struct Node {
    //! @NodeIDStub
    Index id;
    //! Voronoi area associated with the node.
    /**
     * The Voronoi area is the sum of (mixed) Voronoi areas inside the triangles, incident to the node.
     * Definition follows [Gueguen et al. 2017](https://doi.org/10.1039/C7SM01272A).
     * \f[
     * A_{i} = \sum_{j} A'_{ij}.
     * \f]
     *
     * @see Triangulation::mixed_area
     * See Figure tr1. C in Triangulation.
     * @see Node::curvature_vec Triangulation::update_bulk_node_geometry(Index)
     */
    Real area;
    //! If the node is part of a closed surface triangulation, then the `volume` contains the volume of the tetrahedron
    //! connected to each voronoi cell sub-triangle and the center of the lab coordinate system  as defined in [Gueguen
    //! et al. 2017](https://doi.org/10.1039/C7SM01272A).
    /** * This means that the volume of an individual node does not have a proper physical interpretation.
     * Only the sum of all node volumes, which is given by the triangulation
     * is interpretable as a physical volume of an object.
     * The definition follows [Gueguen et al. 2017](https://doi.org/10.1039/C7SM01272A).
     * \f[
     * V_{ij} = A_{ij} \vec{x}_{i}\cdot \frac{\vec{n}_{ij,j+1}}{\| \vec{n}_{ij,j+1} \|}.
     * \f]
     * See Figure tr1. D in Triangulation.
     * @see Node::curvature_vec Triangulation::update_bulk_node_geometry(Index)
     */
    Real volume;
    //! `unit_bending_energy` corresponds to the [Helfrich bending
    //! energy](https://en.wikipedia.org/wiki/Elasticity_of_cell_membranes) with bending rigidity 1 and gaussian bending
    //! stiffness 0.
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
    Real unit_bending_energy;
    //! Position of the node in the lab frame.
    vec3<Real> pos;
    //! Curvature vector of the node.
    /**
     * The definition of the curvature vector follows [Meyer et al. 2003](https://doi.org/10.1007/978-3-662-05105-4_2).
     * \f[
     * \vec{K}_i = \frac{1}{2A_i}\sum_{j(i)} \left( \cot\left(\alpha_{ij}^{j+1}\right) +
     * \cot\left(\alpha_{ij}^{j-1}\right) \right)\vec{\ell}_{ij}
     * \f]
     * See Figure tr1. B in Triangulation.
     * @see Node::curvature_vec Triangulation::update_bulk_node_geometry(Index)
     */
    vec3<Real> curvature_vec;
    //! A vector containing the global ids of the current node's next neighbors.
    /**
     * `nn_ids` contains the ids of nodes that are connected to this node in the triangulation.
     * The next neighbors that are also mutual neighbors in the triangulation are stored sequentially in the vector.
     * The last and the first elements are also neighbors, i.e., the nn_ids vector wraps around.
     * During the calculation, this is facilitated through the use of @ref fp::Neighbors.
     * @note The order of the next neighbors matters for the proper function of fp::Triangulation but is not guaranteed
     * by this data structure. See Figure tr1. A, in Triangulation.
     */
    NeighborsArray<Index> nn_ids;
    //! Distance vectors pointing from the node to its next neighbors.
    NeighborsArray<vec3<Real>> nn_distances;
    //! The Verlet list contains the ids of nodes that are close to this node.
    std::vector<Index> verlet_list;

    // unit-tested
    //! Find and deletes the element with the id `to_pop_nn_id` in the `nn_id` vector.
    void pop_nn(Index to_pop_nn_id) {
        /**
         * @param to_pop_nn_id @NNIDStub This id is supposed to be removed from the next neighbor id vector.
         * @see Node::nn_ids
         * @note this will lead to resizing of the vector, which can be expensive!
         * @warning If the provided next neighbor id is not part of the Node::nn_ids, this function will fail silently.
         * It will not delete anything and won't throw any errors or warnings;
         */
        auto pop_pos = find_nns_loc_pointer(to_pop_nn_id);
        auto dist    = pop_pos - nn_ids.begin();

        if (pop_pos != nn_ids.end()) {
            // I checked that this would work on example code on cppreference https://godbolt.org/z/6qf8c9nTz
            nn_ids.erase(pop_pos);
            nn_distances.erase(nn_distances.begin() + dist);
        }
    }

    auto find_nns_loc_pointer(Index nn_id) {
        /**
         * @brief Given the global id of the next neighbor, this function can be used to locate it in the Node::nn_ids
         * vector.
         *
         * This function is just a convenient wrapper around the
         * [std::find](https://en.cppreference.com/w/cpp/algorithm/find) function.
         * ```
         * std::find(nn_ids.begin(), nn_ids.end(), to_pop_nn_id);
         * ```
         * @param nn_id @NNIDStub
         * @return if `nn_id` is contained in Node::nn_ids then the pointer to the position of that id in the `nn_ids`
         * vector will be returned. Otherwise `nn_ids.end()`.
         * @warning This function is not responsible for graceful handling of `nn_id`'s that are not found in the
         * Node::nn_ids vector. If the `nn_id` is not contained in Node::nn_ids then the `nn_ids.end()` iterator will be
         * returned. It is up to the user to perform the necessary checks to avoid undefined behavior that might result
         * from trying to delete uninitiated memory.
         */
        return std::find(nn_ids.begin(), nn_ids.end(), nn_id);
    }

    // unit-tested
    void emplace_nn_id(Index to_emplace_nn_id, const vec3<Real> &to_emplace_nn_pos, Index loc_idx) {
        /**
         * @brief This function can be used to add new next neighbors to this node.
         *
         * This function constructs `to_emplace_nn_id` right before `to_emplace_pos`,
         * i.e. if to_emplace_nn_id is 3, to_emplace_nn_id will be constructed right before the
         * 3rd element and will become the new 3rd element.
         * @param to_emplace_nn_id @NNIDStub This id is supposed to be added to the Node::nn_ids vector of this node.
         * @param to_emplace_nn_pos const reference to the 3 dimensional position vector (type vec3<Real>) containing
         * the position of the new next neighbour. This input is used to calculate the correct distance between this
         * node and the new next neighbor, which will then be added to the Node::nn_distances vector.
         * @param loc_idx @LocNNIndexStub
         * @note This function causes the resizing of two vectors, which can be costly.
         * @warning Making next neighbors is a symmetric operation. I.e., if node one becomes the next neighbor of node
         * two, node two also has to become the next neighbor of node one. However, this function is not responsible for
         * this relationship. It only adds a new next neighbor to this node, and the higher-order structures, like
         * Triangulation, are responsible for guaranteeing the reciprocal relationship.
         * @see Triangulation::emplace_before(Index, Index, Index)
         */
        if (loc_idx < nn_ids.size()) {
            nn_ids.emplace(loc_idx, to_emplace_nn_id);
            nn_distances.emplace(loc_idx, to_emplace_nn_pos - pos);
        }
    }

    // unit-tested
    //! This function can provide the stored distance vector to the next neighbor.
    const vec3<Real> &get_distance_vector_to(Index nn_id) const {
        /**
         * @param nn_id @NNIDStub.
         * @return returns the distance currently stored in the Node::nn_distances vector for the requested next
         * neighbor. If the provided `nn_id` can not be found in the Node::nn_ids vector, then the function writes an
         * error message to standard error output and terminates the program with exit code 12.
         * @note @TerminationNoteStub
         */
        auto id_pos = std::find(nn_ids.begin(), nn_ids.end(), nn_id);
        if (id_pos != nn_ids.end()) {
            return nn_distances[static_cast<Index>(id_pos - nn_ids.begin())];
        } else {
            std::cerr << "nn_id:" << nn_id
                      << " provided to `get_distance_vector_to` is not a next neighbour of the node " << id;
            exit(12);
        }
    }

    // defaulted operators are not explicitly unit-tested
    /**
     * @brief Default equality operator.
     *
     * @param other_node constant reference to the other Node.
     * @return True if both nodes are equal.
     */
    bool operator==(const Node<Real, Index> &other_node) const = default;

    /**
     * @brief Streaming operator that can print formatted output to standard out with all Node data fields.
     *
     * @param os This is intended to be std::cout or any other std::ofstream reference.
     * @param node The streamed node.
     * @return Updated stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const Node<Real, Index> &node) {

        os << "node: " << node.id << '\n'
           << "area: " << node.area << '\n'
           << "volume: " << node.volume << '\n'
           << "unit_bending_energy: " << node.unit_bending_energy << '\n'
           << "curvature_vec: " << node.curvature_vec << '\n'
           << "pos: " << node.pos << '\n'
           << "nn_ids: ";
        for (const auto &nn_id : node.nn_ids) {
            os << nn_id << ' ';
        }
        os << '\n' << "nn_distances: ";
        for (const auto &nn_dist : node.nn_distances) {
            os << nn_dist << '\n';
        }
        os << '\n';

        return os;
    }
};

} // namespace fp::experimental
#endif // FLIPPY_NODES_HPP
