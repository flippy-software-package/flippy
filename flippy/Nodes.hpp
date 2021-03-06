#ifndef FLIPPY_NODES_HPP
#define FLIPPY_NODES_HPP

#include <vector>
#include <unordered_set>

#include "external/json.hpp"

#include "vec3.hpp"

namespace fp {
using Json = nlohmann::json;

/**
   * This is a DUMB DATA STRUCTURE, meaning that it is not responsible for the coherence of the data that it contains.
   * I.e. it will never check if the curvature is the norm of the curvature vector or if the nn_ids and nn_distances are in the correct order.
   * However it does check the data for consistency. It will match the length of nn_ids and nn_distances. And pop and add both of them together.
   */
template<floating_point_number Real, integer_number Index>
struct Node
{
  //! a number between 0 and max number of nodes - 1
  Index id;
  //! voronoi area associated to the node
  Real area;
  /**
   * if the node is part of a closed surface triangulation, then the `volume` gives the volume
   * of the tetrahedron connected to each voronoi cell sub-triangle and the center of the lab coordinate system  as defined in [1].
   * This means that the volume of an individual node does not have a proper physical interpretation. Only the sum of all node
   * volumes, which is given by the triangulation is interpretable as a physical volume of an object.
   */
  Real volume;
  /**
   * unit_bending energy corresponds to the Helfrich bending energy with bending rigidity 1
   * \f[
   *  \mathrm{unit\_bending\_energy} = \frac{1}{2} A_{\mathrm{node}} (2 H_{node})^2
   * \f]
   * where \f$ H_{node} \f$ is the mean curvature of the node given by:
   * \f[
   * H_{node}^2 = \frac{\vec{K}_{node}}{2A_{node}} \cdot \frac{\vec{K}_{node}}{2A_{node}}
   * \f],
   * with  \f$ \vec{K} \f$ denoting the `curvature_vector` member of the node
   */
  Real unit_bending_energy;
  //! position of the node in the lab frame
  vec3<Real> pos;
  //! curvature vector of the node as defined in [1]
  vec3<Real> curvature_vec;
  //! nn_ids contains the ids of nodes that are connected to this node
  std::vector<Index> nn_ids;
  //! distance vectors pointing from the node to it's next neighbours.
  std::vector<vec3<Real>> nn_distances;
  //! Verlet list contains the ids of nodes that are close to this node
  std::vector<Index> verlet_list;

  // unit-tested
  void pop_nn(Index to_pop_nn_id)
  {
      // finds the element with the id to_pop_nn_id in the nn_id vector and deletes it.
      // this will lead to resizing of the vector!
      auto pop_pos = find_nns_loc_pointer(to_pop_nn_id);
      Index dist = pop_pos - nn_ids.begin();

      if (pop_pos!=nn_ids.end()) {
          // I checked that this would work on example code on cppreference https://godbolt.org/z/6qf8c9nTz
          nn_ids.erase(pop_pos);
          nn_distances.erase(nn_distances.begin() + dist);
      }
  }

  auto find_nns_loc_pointer(Index to_pop_nn_id){
      return std::find(nn_ids.begin(), nn_ids.end(), to_pop_nn_id);
  }
  // unit-tested
  void emplace_nn_id(Index to_emplace_nn_id, vec3<Real> const& to_emplace_nn_pos, Index loc_idx)
  {
      /** this constructs to_emplace_nn_id right before to_emplace_pos
       * ie if to_emplace_nn_id is 3, to_emplace_nn_id will be constructed right before the
       * 3rd element and will become the new 3rd element.
       */
      if (loc_idx<(Index) nn_ids.size()) {
          nn_ids.emplace(nn_ids.begin() + loc_idx, to_emplace_nn_id);
          nn_distances.emplace(nn_distances.begin() + loc_idx, to_emplace_nn_pos - pos);
      }
  }

  //unit-tested
  vec3<Real> const& get_distance_vector_to(Index nn_id) const
  {
      auto id_pos = std::find(nn_ids.begin(), nn_ids.end(), nn_id);
      if (id_pos!=nn_ids.end()) {
          return nn_distances[(id_pos - nn_ids.begin())];
      }
      else {
          std::cerr << "nn_id:" << nn_id << " provided to `get_distance_vector_to` is not a next neighbour of the node "
                    << id;
          exit(12);
      }
  }

  bool operator==(Node<Real, Index> const& other_node) const = default;

  friend std::ostream& operator<<(std::ostream& os, Node<Real, Index> const& node1)
  {
      os << "node: " << node1.id << '\n'
         << "area: " << node1.area << '\n'
         << "volume: " << node1.volume << '\n'
         << "unit_bending_energy: " << node1.unit_bending_energy << '\n'
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

template<floating_point_number Real, integer_number Index>
struct Nodes
{
    std::vector<Node<Real, Index>> data;
    Nodes() = default;
    explicit Nodes(std::vector<Node<Real, Index>> data_inp)
            :data(data_inp) { }
    explicit Nodes(Json const& node_dict)
    {
        /**
         * Initiating nodes from a Json of a node collection.
         * The nodes in the json file must be sequentially numbered from 0 to Number_of_nodes - 1
         */
        std::vector<Index> nn_ids_temp, verlet_list_temp;
        data.resize((node_dict.size()));
        for (auto const& node: node_dict.items()) {
            auto const& node_id = node.key();
            Index node_index = std::stol(node_id);
            auto const& raw_pos = node.value()["pos"];
            vec3<Real> pos{(Real) raw_pos[0], (Real) raw_pos[1], (Real) raw_pos[2]};

            auto const& raw_curv = node.value()["curvature_vec"];
            vec3<Real> curvature_vec{(Real) raw_curv[0], (Real) raw_curv[1], (Real) raw_curv[2]};
            Real unit_bending_energy = node.value()["unit_bending_energy"];
            Real area = node.value()["area"];
            Real volume = node.value()["volume"];

            nn_ids_temp = node_dict[node_id]["nn_ids"].get<std::vector<Index>>();
            verlet_list_temp = node_dict[node_id]["verlet_list"].get<std::vector<Index>>();
            std::vector<vec3<Real>>
                    nn_distances;

            data[node_index] = Node<Real, Index>{
                    .id{node_index},
                    .area{area},
                    .volume{volume},
                    .unit_bending_energy{unit_bending_energy},
                    .pos{pos},
                    .curvature_vec{curvature_vec},
                    .nn_ids{nn_ids_temp},
                    .nn_distances{nn_distances},
                    .verlet_list{verlet_list_temp}
            };
        }
    }

    typename std::vector<Node<Real, Index>>::iterator begin(){return data.begin();}
    typename std::vector<Node<Real, Index>>::const_iterator begin() const {return data.begin();}
    typename std::vector<Node<Real, Index>>::iterator end() {return data.end();}
    typename std::vector<Node<Real, Index>>::const_iterator end() const {return data.end();}

    // getters and setters
    //unit-tested
    const vec3<Real>& pos(Index node_id) const {return data[node_id].pos;}
    //unit-tested
    void set_pos(Index node_id, vec3<Real> const& new_pos){data[node_id].pos=new_pos;}
    void set_pos(Index node_id, vec3<Real> && new_pos){data[node_id].pos=new_pos;}
    void displace(Index node_id, vec3<Real>const& displ){data[node_id].pos+=displ;}
    void displace(Index node_id, vec3<Real>&& displ){data[node_id].pos+=displ;}

    const vec3<Real>& curvature_vec(Index node_id) const {return data[node_id].curvature_vec;}
    void set_curvature_vec(Index node_id, vec3<Real> const& new_cv) {data[node_id].curvature_vec=new_cv;}
    void set_curvature_vec(Index node_id, vec3<Real> && new_cv) {data[node_id].curvature_vec=new_cv;}

    Real area(Index node_id)const{return data[node_id].area;}
    void set_area(Index node_id, Real new_area){data[node_id].area = new_area;}

    Real volume(Index node_id)const{return data[node_id].volume;}
    void set_volume(Index node_id, Real new_volume){data[node_id].volume = new_volume;}

    Real unit_bending_energy(Index node_id)const{return data[node_id].unit_bending_energy;}
    void set_unit_bending_energy(Index node_id, Real new_sce){data[node_id].unit_bending_energy=new_sce;}

    //unit-tested
    const auto& nn_ids(Index node_id)const{return data[node_id].nn_ids;}
    //unit-tested
    void set_nn_ids(Index node_id, std::vector<Index>const& new_nn_ids){data[node_id].nn_ids = new_nn_ids;}
    //unit-tested
    Index nn_id(Index node_id, Index loc_nn_index)const{return data[node_id].nn_ids[loc_nn_index];}
    //unit-tested
    void set_nn_id(Index node_id, Index loc_nn_index, Index nn_id){data[node_id].nn_ids[loc_nn_index]=nn_id;}
    void emplace_nn_id(Index node_id, Index to_emplace_nn_id, Index loc_idx)
    {data[node_id].emplace_nn_id(to_emplace_nn_id, pos(to_emplace_nn_id), loc_idx);}

    const auto& nn_distances(Index node_id)const{return data[node_id].nn_distances;}
    const auto& get_nn_distance_vector_between(Index node_id, Index nn_id) const{
        return data[node_id].get_distance_vector_to(nn_id);
    }
    void set_nn_distance(Index node_id, Index loc_nn_index, vec3<Real>&& dist){data[node_id].nn_distances[loc_nn_index]=dist;}
    void set_nn_distance(Index node_id, Index loc_nn_index, vec3<Real> const& dist){data[node_id].nn_distances[loc_nn_index]=dist;}

    [[nodiscard]] Index size() const { return data.size(); }
    Node<Real, Index>& operator[](Index idx) { return data[idx]; }
    const Node<Real, Index>& operator[](Index idx) const { return data.at(idx); }

    [[nodiscard]] Json make_data() const
    {
        Json json_data;
        for (auto& node : data) {
            json_data[std::to_string(node.id)] = {
                    {"area", node.area},
                    {"volume", node.volume},
                    {"unit_bending_energy", node.unit_bending_energy},
                    {"pos", {node.pos[0], node.pos[1], node.pos[2]}},
                    {"curvature_vec", {node.curvature_vec[0], node.curvature_vec[1], node.curvature_vec[2]}},
                    {"nn_ids", node.nn_ids},
                    {"verlet_list", node.verlet_list},
            };
        }
        return json_data;
    }
};
}
#endif //FLIPPY_NODES_HPP
