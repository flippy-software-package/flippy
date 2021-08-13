#ifndef EXPERIMENT_H__NODES_H_
#define EXPERIMENT_H__NODES_H_

#include <vector>
#include <unordered_set>

#include "external/json.hpp"

#include "vec3.hpp"
#include "utils.h"

namespace fp {
using Json = nlohmann::json;

template<typename T>
static bool is_member(std::vector<T> const& v, T const& el)
{
    /***
     * if the function returns true than el is contained in v (at least once).
     */
    return (std::find(v.begin(), v.end(), el)!=v.end());
}

template<typename Real, typename Index>
struct Node
{
  /**
   * This is a DUMB DATA STRUCTURE, meaning that it is not responsible for the coherence of the data that it contains.
   * I.e. it will never check if the curvature is the norm of the curvature vector or if the nn_ids and nn_distances are in the correct order.
   * However it does check the data for consistency. It will match the length of nn_ids and nn_distances. And pop and add both of them together. //Todo this is not yet implemented
   */
  Index id;
  Real area;
  Real volume;
  Real scaled_curvature_energy;
  vec3<Real> pos;
  vec3<Real> curvature_vec;
  std::vector<Index> nn_ids;
  std::vector<Index> verlet_list;
  std::vector<vec3<Real>> nn_distances;

  //Todo  unittest
  void pop_verlet_neighbor(Index const& to_pop_vn_id)
  {
      // finds the element with the id to_pop_nn_id in the nn_id vector and deletes it.
      // this will lead to resizing of the vector!
      auto pop_pos = std::find(nn_ids.begin(), nn_ids.end(), to_pop_vn_id);

      if (pop_pos!=nn_ids.end()) {
          // I checked that this would work on example code on cppreference https://godbolt.org/z/6qf8c9nTz
          verlet_list.erase(pop_pos);
      }
  }

  // unit-tested
  void pop_nn(Index const& to_pop_nn_id)
  {
      // finds the element with the id to_pop_nn_id in the nn_id vector and deletes it.
      // this will lead to resizing of the vector!
      auto pop_pos = std::find(nn_ids.begin(), nn_ids.end(), to_pop_nn_id);
      Index dist = pop_pos - nn_ids.begin();

      if (pop_pos!=nn_ids.end()) {
          // I checked that this would work on example code on cppreference https://godbolt.org/z/6qf8c9nTz
          nn_ids.erase(pop_pos);
          nn_distances.erase(nn_distances.begin() + dist);
      }
  }

  // unit-tested
  void emplace_nn_id(Index const& to_emplace_nn_id, vec3<Real> const& to_emplace_nn_pos, Index const& to_emplace_pos)
  {
      /** this constructs to_emplace_nn_id right before to_emplace_pos
       * ie if to_emplace_nn_id is 3, to_emplace_nn_id will be constructed right before the
       * 3rd element and will become the new 3rd element.
       */
      if (to_emplace_pos<(Index) nn_ids.size()) {
          nn_ids.emplace(nn_ids.begin() + to_emplace_pos, to_emplace_nn_id);
          nn_distances.emplace(nn_distances.begin() + to_emplace_pos, to_emplace_nn_pos - pos);
      }
  }

  vec3<Real> const& get_distance_vector_to(const Index& nn_id) const
  {
      auto id_pos = std::find(nn_ids.begin(), nn_ids.end(), nn_id);
      if (id_pos!=nn_ids.end()) {
          return nn_distances[(Index) (id_pos - nn_ids.begin())];
      }
      else {
          std::cerr << "nn_id:" << nn_id << " provided to `get_distance_vector_to` is not a next neighbour of the node"
                    << id;
          exit(12);
      }
  }

  void copy_only_geometry(Node<Real, Index> const& node_to_copy)
  {
      area = node_to_copy.area;
      volume = node_to_copy.volume;
      scaled_curvature_energy = node_to_copy.scaled_curvature_energy;
      pos = node_to_copy.pos;
      curvature_vec = node_to_copy.curvature_vec;
      nn_distances = node_to_copy.nn_distances;
  }

  bool operator==(Node<Real, Index> const& other_node) const = default;

  friend std::ostream& operator<<(std::ostream& os, Node<Real, Index> const& node1)
  {
      os << "node: " << node1.id << '\n'
         << "area: " << node1.area << '\n'
         << "volume: " << node1.volume << '\n'
         << "scaled_curvature_energy: " << node1.scaled_curvature_energy << '\n'
         << "curvature_vec: " << node1.curvature_vec << '\n'
         << "pos: " << node1.pos << '\n'
         << "nn_ids: ";
      for (auto const& nn_id: node1.nn_ids) {
          os << nn_id << ' ';
      }
      os << '\n'
         << "nn_distances: ";
      for (auto const& nn_dist: node1.nn_distances) {
          os << nn_dist << '\n';
      }
      os << '\n';

      return os;
  }

};

template<typename Real, typename Index>
struct Nodes
{
    std::vector<Node<Real, Index>> data;
    Nodes() = default;
    Nodes(std::vector<Node<Real, Index>> data_inp, Real const& verlet_radius_inp)
            :data(data_inp), verlet_radius(verlet_radius_inp), verlet_radius_squared(verlet_radius*verlet_radius) { }
    explicit Nodes(Json const& node_dict, Real const& verlet_radius_inp)
            :verlet_radius(verlet_radius_inp)
    {
        /*
         * Initiating nodes from a Json of a node collection.
         * The nodes in the json file must be sequentially numbered from 0 to Number_of_nodes - 1
         */
        verlet_radius_squared = verlet_radius*verlet_radius;
        std::vector<Index> nn_ids_temp;
        data.resize((node_dict.size()));
        for (auto const& node: node_dict.items()) {
            auto const& node_id = node.key();
            Index node_index = std::stol(node_id);
            auto const& raw_pos = node.value()["pos"];
            vec3<Real> pos{(Real) raw_pos[0], (Real) raw_pos[1], (Real) raw_pos[2]};

            auto const& raw_curv = node.value()["curvature_vec"];
            vec3<Real> curvature_vec{(Real) raw_curv[0], (Real) raw_curv[1], (Real) raw_curv[2]};
            Real scaled_curvature_energy = node.value()["scaled_curvature_energy"];
            Real area = node.value()["area"];
            Real volume = node.value()["volume"];

            nn_ids_temp = std::vector<Index>(node_dict[node_id]["nn_ids"]);
            std::vector<vec3<Real>>
            nn_distances;

            data[node_index] = (Node<Real, Index>{
                    .id{node_index},
                    .area{area},
                    .volume{volume},
                    .scaled_curvature_energy{scaled_curvature_energy},
                    .pos{pos},
                    .curvature_vec{curvature_vec},
                    .nn_ids{nn_ids_temp},
//                    .verlet_list={std::vector<Index>(10)},
                    .nn_distances{nn_distances}});
        }
    }

    void make_verlet_list()
    {
        for (auto& node: data) {
            node.verlet_list.clear();
        }
        for (auto node_p = data.begin(); node_p!=data.end(); ++node_p) {
            for (auto other_node_p = data.begin(); other_node_p!=node_p; ++other_node_p) {
                if ((node_p->pos - other_node_p->pos).norm_square()<verlet_radius_squared) {
                    node_p->verlet_list.push_back(other_node_p->id);
                    other_node_p->verlet_list.push_back(node_p->id);
                }
            }

        }

    }

    [[nodiscard]] Index size() const { return data.size(); }
    Node<Real, Index>& operator[](Index const& idx) { return data[idx]; }
    const Node<Real, Index>& operator[](Index const& idx) const { return data.at(idx); }

    [[nodiscard]] Json make_data() const
    {
        Json json_data;
        for (auto& node : data) {
            json_data[std::to_string(node.id)] = {
                    {"area", node.area},
                    {"volume", node.volume},
                    {"scaled_curvature_energy", node.scaled_curvature_energy},
                    {"pos", {node.pos[0], node.pos[1], node.pos[2]}},
                    {"curvature_vec", {node.curvature_vec[0], node.curvature_vec[1], node.curvature_vec[2]}},
                    {"nn_ids", node.nn_ids},
            };
        }
        return json_data;
    }
private:
    Real verlet_radius;
    Real verlet_radius_squared;
};
}
#endif //EXPERIMENT_H__NODES_H_
