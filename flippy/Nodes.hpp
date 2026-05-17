/*
 *```txt
 *
 *  .d888 888 d8b
 * d88P"  888 Y8P
 * 888    888
 * 888888 888 888 88888b.  88888b.  888  888
 * 888    888 888 888 "88b 888 "88b 888  888     simulating package for
 * 888    888 888 888  888 888  888 888  888     dynamically triangulated
 * 888    888 888 888 d88P 888 d88P Y88b 888     surfaces
 * 888    888 888 88888P"  88888P"   "Y88888
 *                888      888           888
 *                888      888      Y8b d88P
 *                888      888       "Y88P"
 *
 * https://github.com/flippy-software-package/flippy
 *
 *
 * MIT License
 *
 * Copyright (c) 2021 George Dadunashvili
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *```
 */

#ifndef FLIPPY_NODES_HPP
#define FLIPPY_NODES_HPP
/**
 * @file
 * @brief This file contains the fp::Node and fp::Nodes classes, data structures that represent a single node of the
 * triangulation and the collection of all nodes of the triangulation, respectively.
 */
#include "external/json.hpp"
#include "vec3.hpp"

#include <unordered_set>
#include <utility>
#include <vector>

namespace fp {
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
struct Node {
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
    /**
     * This means that the volume of an individual node does not have a proper physical interpretation.
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
    std::vector<Index> nn_ids;
    //! Distance vectors pointing from the node to its next neighbors.
    std::vector<vec3<Real>> nn_distances;
    //! The Verlet list contains the ids of nodes that are close to this node.
    std::vector<Index> verlet_list;

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
            auto signed_loc_idx = static_cast<long long>(loc_idx);
            nn_ids.emplace(nn_ids.begin() + signed_loc_idx, to_emplace_nn_id);
            nn_distances.emplace(nn_distances.begin() + signed_loc_idx, to_emplace_nn_pos - pos);
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
    bool operator==(const Node &other_node) const = default;

    /**
     * @brief Streaming operator that can print formatted output to standard out with all Node data fields.
     *
     * @param os This is intended to be std::cout or any other std::ofstream reference.
     * @param node The streamed node.
     * @return Updated stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const Node &node) {

        os << "node: " << node.id << '\n'
           << "area: " << node.area << '\n'
           << "volume: " << node.volume << '\n'
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

/**
 * @brief Data structure containing all nodes of the Triangulation.
 *
 * The Nodes struct is capable of reinitializing nodes from a well-formed JSON object or from a simple
 * [std::vector](https://en.cppreference.com/w/cpp/container/vector) that contains all nodes of a triangulation. The
 * nodes class is basically a wrapper around a vector of nodes, i.e., `std::vector<Node>`, and provides additional
 * functionality to manipulate and query this data structure. Nodes class is also meant to be the interface with which
 * the end user is manipulating individual nodes.
 * @tparam Real @RealStub
 * @tparam Index @IndexStub
 */
struct Nodes {
    std::vector<Node> data; //!< Data member that contains the individual nodes.

    Nodes() = default; //!< Default constructor.
    explicit Nodes(std::vector<Node> data_inp) : data(data_inp) {
        /**
         * Copies the data from a vector of nodes and creates a new Nodes struct.
         * @param data_inp A standard vector containing all the nodes that are supposed to create a new Nodes class.
         */
    } //!< Constructor from a vector.
    explicit Nodes(const Json &node_dict) {
        /**
         * Initiating nodes from a JSON object of a node collection.
         * The nodes in the JSON file must be sequentially numbered from 0 to Number_of_nodes - 1.
         * @param node_dict JSON object that contains a collection of nodes.
         * @warning If the JSON object is malformed, then the constructor will fail and propagate a runtime error from
         * the JSON parser.
         */
        std::vector<Index> nn_ids_temp, verlet_list_temp;
        data.resize((node_dict.size()));
        for (const auto &node : node_dict.items()) {
            const auto             &node_id    = node.key();
            fp::IndexingNumber auto node_index = static_cast<Index>(std::stol(node_id));
            const auto             &raw_pos    = node.value()["pos"];
            vec3<Real>              pos{(Real)raw_pos[0], (Real)raw_pos[1], (Real)raw_pos[2]};

            const auto &raw_curv = node.value()["curvature_vec"];
            vec3<Real>  curvature_vec{(Real)raw_curv[0], (Real)raw_curv[1], (Real)raw_curv[2]};
            Real        area   = node.value()["area"];
            Real        volume = node.value()["volume"];

            nn_ids_temp      = node_dict[node_id]["nn_ids"].get<std::vector<Index>>();
            verlet_list_temp = node_dict[node_id]["verlet_list"].get<std::vector<Index>>();
            std::vector<vec3<Real>> nn_distances;

            data[static_cast<size_t>(node_index)] = Node{.id     = node_index,
                                                         .area   = area,
                                                         .volume = volume,
                                                         .pos{pos},
                                                         .curvature_vec{curvature_vec},
                                                         .nn_ids{nn_ids_temp},
                                                         .nn_distances{nn_distances},
                                                         .verlet_list{verlet_list_temp}};
        }
    } //!< Constructor from JSON.

    typename std::vector<Node>::iterator begin() {
        /**
         * This function allows the Nodes struct to be used in range-based `for` loops.
         * @return `data.begin()`
         */
        return data.begin();
    } //!< Returns an iterator to the beginning of the underlying data member that contains the collection of the
      //!< nodes.
    [[nodiscard]] typename std::vector<Node>::const_iterator begin() const {
        /**
         * This function allows the Nodes struct to be used in range-based `for` loops in constant environments.
         * @return a constant iterator `data.begin()`.
         */
        return data.begin();
    } //!< \overload

    typename std::vector<Node>::iterator end() {
        /**
         * This function allows the Nodes struct to be used in range-based `for` loops.
         * @return `data.end()`.
         */
        return data.end();
    } //!< Returns an iterator to the end of the underlying data member that contains the collection of the nodes.
    [[nodiscard]] typename std::vector<Node>::const_iterator end() const {
        /**
         * This function allows the Nodes struct to be used in range-based `for` loops in constant environments.
         * @return a constant iterator `data.end()`.
         */
        return data.end();
    } //!< \overload

    // getters and setters

    // unit-tested

    void emplace_nn_id(Index node_id, Index to_emplace_nn_id, Index loc_nn_index) {
        /**
         * This function is a wrapper around Node::emplace_nn_id(Index , vec3<Real> const& , Index).
         * @param node_id @NodeIDStub
         * @param to_emplace_nn_id @NNIDStub
         * @param loc_nn_index @LocNNIndexStub
         */
        data[node_id].emplace_nn_id(to_emplace_nn_id, data[to_emplace_nn_id].pos, loc_nn_index);
    } //!< Emplace a the id of a new node in the Node::nn_ids vector, in front of the loc_idx position.

    void set_nn_distance(Index node_id, Index loc_nn_index, vec3<Real> &&dist) {
        /**
         * @param node_id @NNIDStub
         * @param loc_nn_index @LocNNIndexStub
         * @param dist rvalue reference to a 3D distance vector (that points from node_id to its next neighbour).
         */
        data[node_id].nn_distances[loc_nn_index] = dist;
    } //!< Overwrite the next neighbor distance with a new 3d vector.

    void set_nn_distance(Index node_id, Index loc_nn_index, const vec3<Real> &dist) {
        /**
         * @param node_id @NNIDStub
         * @param loc_nn_index @LocNNIndexStub
         * @param dist lvalue constant reference to a 3D distance vector (that points from node_id to its next
         * neighbour).
         */
        data[node_id].nn_distances[loc_nn_index] = dist;
    } //!< \overload

    [[nodiscard]] Index size() const {
        return static_cast<Index>(data.size());
    } //!< Size of the Nodes data member. @return Size of the data vector, same as the number of nodes.

    Node &operator[](Index node_id) {
        /**
         * Nodes[node_id] is the same as Nodes.data[node_id].
         * @param node_id @NodeIDStub
         * @return Reference to the Node struct with the id corresponding to node_id.
         */
        return data[node_id];
    } //!< Square bracket operator overload for convenient indexing of the Nodes struct.
    const Node &operator[](Index node_id) const {
        /**
         * Nodes[node_id] in the constant environment is the same as Nodes.data.at(node_id).
         * @param node_id @NodeIDStub
         * @return Constant reference to the Node struct with the id corresponding to node_id.
         */
        return data.at(node_id);
    } //!< @overload

    [[nodiscard]] Json make_data() const {
        /**
         * @return JSON object that represents a serialization of the data contained in Nodes.
         * This JSON object can later be used to reconstruct the Nodes object.
         */
        Json json_data;
        for (auto &node : data) {
            json_data[std::to_string(node.id)] = {
                {"area", node.area},
                {"volume", node.volume},
                {"pos", {node.pos[0], node.pos[1], node.pos[2]}},
                {"curvature_vec", {node.curvature_vec[0], node.curvature_vec[1], node.curvature_vec[2]}},
                {"nn_ids", node.nn_ids},
                {"verlet_list", node.verlet_list},
            };
        }
        return json_data;
    } //!< Serialize the Nodes struct to a JSON object.
};
} // namespace fp
#endif // FLIPPY_NODES_HPP
