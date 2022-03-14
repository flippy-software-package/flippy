#ifndef FLIPPY_TRIANGULATOR_HPP
#define FLIPPY_TRIANGULATOR_HPP

#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "custom_concepts.hpp"
#include "vec3.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846	/* pi */
#endif

/**
 * The API stability of the functions in the implementation namespace is not guaranteed!
 * Functions that are part of the implementation namespace are not part of the public facing API and are not intended fot the end-user.
 * Since flippy is a headers only library this could not be hidden in source files.
 */
namespace fp::implementation{
template<floating_point_number Real, integer_number Index>
struct SimpleNodeData{
  std::string hash{};
  Index id{};
  vec3<Real> pos{};
  std::unordered_set<std::string> nn_hashes{};
};

template<floating_point_number Real, integer_number Index>
class IcosahedronSubTriangulation
{
public:
    static std::string hash_node(Index c)
    {
        /**
         * returns a unique hash for a corner node, which is just the id of that corner node.
         */
        return std::to_string(c);
    }

    static std::string hash_node(Index c0, Index c1, Index n)
    {
        /**
         * returns a unique hash for a node on one of the sides of the initial triangle.
         * This hash is determined by the (ordered) corner nodes of the initial edge and the index of the node.
         */
        Index a = std::min(c0, c1);
        Index b = std::max(c0, c1);
        return std::to_string(a) + "_" + std::to_string(b) + "_" + std::to_string(n);
    }

    static std::string hash_node(Index c0, Index c1, Index c2, Index i, Index j)
    {
        /**
         * returns a unique hash for a node in the bulk of the subtriangulation.
         */
        std::vector<Index> cv{c0, c1, c2};
        std::sort(cv.begin(), cv.end());

        return std::to_string(cv[0]) + "_" + std::to_string(cv[1]) + "_" + std::to_string(cv[2])
                + "_" + std::to_string(i) + "_" + std::to_string(j);
    }

    static vec3<Real> r_S1(Real R, Real t, Real f) {
        vec3<Real> r{R * std::sin(t) * std::cos(f), R * std::sin(t) * std::sin(f), R * std::cos(t)};
        return r;
    }

    static constexpr int N_ICOSA_FACEs = 20;
    static constexpr int N_ICOSA_EDGEs = 30;
    static constexpr int N_ICOSA_NODEs = 12;

    static constexpr std::array<int, N_ICOSA_FACEs> FACE_IDs{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    static constexpr std::array<int, N_ICOSA_NODEs> NODE_IDs{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    static constexpr std::array<std::array<int, 3>, N_ICOSA_FACEs> FACE_CORNER_NODES_ = {
            std::array<int, 3>{NODE_IDs[0], NODE_IDs[5], NODE_IDs[1]},
            std::array<int, 3>{NODE_IDs[0], NODE_IDs[1], NODE_IDs[2]},
            std::array<int, 3>{NODE_IDs[0], NODE_IDs[2], NODE_IDs[3]},
            std::array<int, 3>{NODE_IDs[0], NODE_IDs[3], NODE_IDs[4]},
            std::array<int, 3>{NODE_IDs[0], NODE_IDs[4], NODE_IDs[5]},
            std::array<int, 3>{NODE_IDs[10], NODE_IDs[5], NODE_IDs[4]},
            std::array<int, 3>{NODE_IDs[10], NODE_IDs[6], NODE_IDs[5]},
            std::array<int, 3>{NODE_IDs[6], NODE_IDs[1], NODE_IDs[5]},
            std::array<int, 3>{NODE_IDs[6], NODE_IDs[7], NODE_IDs[1]},
            std::array<int, 3>{NODE_IDs[7], NODE_IDs[2], NODE_IDs[1]},
            std::array<int, 3>{NODE_IDs[7], NODE_IDs[8], NODE_IDs[2]},
            std::array<int, 3>{NODE_IDs[8], NODE_IDs[3], NODE_IDs[2]},
            std::array<int, 3>{NODE_IDs[8], NODE_IDs[9], NODE_IDs[3]},
            std::array<int, 3>{NODE_IDs[9], NODE_IDs[4], NODE_IDs[3]},
            std::array<int, 3>{NODE_IDs[9], NODE_IDs[10], NODE_IDs[4]},
            std::array<int, 3>{NODE_IDs[11], NODE_IDs[10], NODE_IDs[9]},
            std::array<int, 3>{NODE_IDs[11], NODE_IDs[6], NODE_IDs[10]},
            std::array<int, 3>{NODE_IDs[11], NODE_IDs[7], NODE_IDs[6]},
            std::array<int, 3>{NODE_IDs[11], NODE_IDs[8], NODE_IDs[7]},
            std::array<int, 3>{NODE_IDs[11], NODE_IDs[9], NODE_IDs[8]}
    };

    static constexpr const std::array<std::array<int, 3>, N_ICOSA_FACEs> FACE_CORNER_NODES = {
            std::array<int, 3>{0, 5, 1},
            std::array<int, 3>{0, 1, 2},
            std::array<int, 3>{0, 2, 3},
            std::array<int, 3>{0, 3, 4},
            std::array<int, 3>{0, 4, 5},
            std::array<int, 3>{10, 5, 4},
            std::array<int, 3>{10, 6, 5},
            std::array<int, 3>{6, 1, 5},
            std::array<int, 3>{6, 7, 1},
            std::array<int, 3>{7, 2, 1},
            std::array<int, 3>{7, 8, 2},
            std::array<int, 3>{8, 3, 2},
            std::array<int, 3>{8, 9, 3},
            std::array<int, 3>{9, 4, 3},
            std::array<int, 3>{9, 10, 4},
            std::array<int, 3>{11, 10, 9},
            std::array<int, 3>{11, 6, 10},
            std::array<int, 3>{11, 7, 6},
            std::array<int, 3>{11, 8, 7},
            std::array<int, 3>{11, 9, 8}
    };

    static std::unordered_map<std::string, SimpleNodeData<Real, Index>> make_corner_nodes()
    {

        Real R = 1.;
        std::unordered_map<std::string, SimpleNodeData<Real, Index>> base_nodes(N_ICOSA_NODEs);
        base_nodes[hash_node(0)] = {.hash=hash_node(0), .id=0, .pos=r_S1(R, 0., 0.)};
        std::string hash;
        hash.reserve(2);
        for (Index i = 1; i<6; ++i) {
            hash = hash_node(i);
            base_nodes[hash] = {
                    .hash=hash,
                    .id=i,
                    .pos=r_S1(R, M_PI/2. - std::atan(0.5), 2*M_PI*(i - 1.)/5.)};
        }

        for (Index i = 6; i<N_ICOSA_NODEs - 1; ++i) {
            hash = hash_node(i);
            base_nodes[hash] = {
                    .hash=hash,
                    .id=i,
                    .pos=r_S1(R, M_PI/2. + std::atan(0.5), 2*M_PI*(i - 6.5)/5.)};
        }
        hash = hash_node(N_ICOSA_NODEs - 1);
        base_nodes[hash] = {
                .hash=hash,
                .id=static_cast<Index>(N_ICOSA_NODEs - 1),
                .pos=r_S1(R, M_PI, 0.)};
        return base_nodes;
    }

    enum TriangleRegion
    {
      TOP_CORNER, BOTTOM_LEFT_CORNER, BOTTOM_RIGHT_CORNER, LEFT_EDGE, BOTTOM_EDGE, DIAGONAL_EDGE, BULK
    };

    static TriangleRegion get_region(Index i, Index j, Index sizeMinOne)
    {
        if (i==0) { return TOP_CORNER; }
        else if (j==0 && i==sizeMinOne) { return BOTTOM_LEFT_CORNER; }
        else if (j==sizeMinOne && i==sizeMinOne) { return BOTTOM_RIGHT_CORNER; }
        else if (j==0) { return LEFT_EDGE; }
        else if (i==sizeMinOne) { return BOTTOM_EDGE; }
        else if (i==j) { return DIAGONAL_EDGE; }
        else { return BULK; }
    }

    static std::string hash_any(Index c0, Index c1, Index c2, Index i, Index j, Index maxIdx)
    {
        switch (get_region(i, j, maxIdx)) {
            case TOP_CORNER:return hash_node(c0);
            case BOTTOM_LEFT_CORNER:return hash_node(c1);
            case BOTTOM_RIGHT_CORNER:return hash_node(c2);
            case LEFT_EDGE:return hash_node(c0, c1, i);
            case BOTTOM_EDGE:return hash_node(c1, c2, j);
            case DIAGONAL_EDGE:return hash_node(c0, c2, j);
            case BULK:return hash_node(c0, c1, c2, i, j);
            default:
                std::cerr<<"something went wrong! provided indices i: "
                         <<i<<" and j: "
                         <<j<<" together with the maxIdx: "<<maxIdx
                         <<" produced a wrong region.\n";
                exit(12);
        }
    }

    static std::vector<std::string> neighbour_hash_vec(Index c0, Index c1, Index c2, Index i, Index j, Index maxIdx)
    {
        std::vector<std::string> neighbour_hash;
        neighbour_hash.reserve(6);
        switch (get_region(i, j, maxIdx)) {

            case TOP_CORNER:neighbour_hash.push_back(hash_any(c0, c1, c2, 1, 0, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, 1, 1, maxIdx));
                return neighbour_hash;

            case BOTTOM_LEFT_CORNER:neighbour_hash.push_back(hash_any(c0, c1, c2, i, j + 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j, maxIdx));
                return neighbour_hash;

            case BOTTOM_RIGHT_CORNER:neighbour_hash.push_back(hash_any(c0, c1, c2, i, j - 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j - 1, maxIdx));
                return neighbour_hash;

            case LEFT_EDGE:neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i + 1, j, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i, j + 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i + 1, j + 1, maxIdx));
                return neighbour_hash;

            case BOTTOM_EDGE:neighbour_hash.push_back(hash_any(c0, c1, c2, i, j - 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j - 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i, j + 1, maxIdx));
                return neighbour_hash;

            case DIAGONAL_EDGE:neighbour_hash.push_back(hash_any(c0, c1, c2, i, j - 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j - 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i + 1, j, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i + 1, j + 1, maxIdx));
                return neighbour_hash;

            case BULK:neighbour_hash.push_back(hash_any(c0, c1, c2, i, j - 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j - 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i - 1, j, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i, j + 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i + 1, j + 1, maxIdx));
                neighbour_hash.push_back(hash_any(c0, c1, c2, i + 1, j, maxIdx));
                return neighbour_hash;
            default:
                std::cerr<<"something went wrong! provided indices i: "
                         <<i<<" and j: "
                         <<j<<" together with the maxIdx: "<<maxIdx
                         <<" produced a wrong region.\n";
                exit(12);
        }

    }

    static Real even_angular_distance_length(Real l, Index k, Index n, Real R = 1.)
    {
        /**
         * The points of the sub-triangulation can not be equally spaced or their angular distances
         * won't be the same.
         */
        if (k==0) {
            return 0;
        }
        else {
            Real fr = static_cast<Real>(k)/static_cast<Real>(n);
            Real denominator = l + sqrt(4.*R*R - l*l)/tan(fr*2.*asin(l/(2.*R)));
            return 2.*R*R/denominator;
        }
    }

    static vec3<Real> get_pos(vec3<Real> const& p0, vec3<Real> const& p1, vec3<Real> const& p2, Index i, Index j, Index maxIdx)
    {
/**
 * get the position o a node in the sub triangulation of a face of the initial icosahedron.
 * ```{.txt}
 *              p0
 *            /___\
 *     e1   /__\/__\ e2
 *        / \ / \  / \
 *      p1 --------- p3
 *            e3
 * ```
 */
        vec3<Real> e1 = p1 - p0;
        vec3<Real> e2 = p2 - p0;
        auto e = e1.norm();
        Real wi = even_angular_distance_length(e, i, maxIdx);

        e1.normalize();
        e2.normalize();

        vec3<Real> li{};
        Real li_norm{0};
        vec3<Real> interm_1 = p0 + wi*e1;
        interm_1.normalize();
        if (i!=0) {
            vec3<Real> interm_2 = p0 + wi*e2;
            interm_2.normalize();

            li = interm_2 - interm_1;
            li_norm = li.norm();
            li.normalize();
        }
        Real wj = even_angular_distance_length(li_norm, j, i);
        return interm_1 + wj*li;
    }

    static std::tuple<Index, Index, Index> get_sorted_face_nodes(std::array<int, 3> face)
    {
        std::sort(face.begin(), face.end());
        return {static_cast<Index>(face[0]),
                static_cast<Index>(face[1]),
                static_cast<Index>(face[2])};
    }
    static void make_face_nodes(std::unordered_map<std::string, SimpleNodeData<Real, Index>>& node_cache, Index nIter)
    {
        vec3<Real> p0, p1, p2, pos;
        Index nEdge = nIter + 2;// total Number of nodes on an edge
        Index maxIdx = nIter + 1;// max value i or j can have

        std::string hash;
        hash.reserve(10);

        std::string c0_h, c1_h, c2_h;
        for (auto face: FACE_CORNER_NODES) {
            auto[c0, c1, c2] = get_sorted_face_nodes(face);
            c0_h = hash_node(c0);
            c1_h = hash_node(c1);
            c2_h = hash_node(c2);
            p0 = node_cache[c0_h].pos;
            p1 = node_cache[c1_h].pos;
            p2 = node_cache[c2_h].pos;
            for (Index i = 0; i<nEdge; ++i) {
                for (Index j = 0; j<=i; ++j) {
                    pos = get_pos(p0, p1, p2, i, j, maxIdx);
                    hash = hash_any(c0, c1, c2, i, j, maxIdx);
                    node_cache[hash] = {.hash=hash, .pos=pos};
                }
            }
        }

        for (auto face: FACE_CORNER_NODES) {
//        std::sort(face.begin(),face.end());
//        auto [c0, c1, c2] = face;
            auto[c0, c1, c2] = get_sorted_face_nodes(face);
            for (Index i = 0; i<nEdge; ++i) {
                for (Index j = 0; j<=i; ++j) {
                    hash = hash_any(c0, c1, c2, i, j, maxIdx);
                    std::vector<std::string> neighbour_hashes = neighbour_hash_vec(c0, c1, c2, i, j, maxIdx);
                    for (auto const& neighbour_hash: neighbour_hashes) {
                        node_cache[hash].nn_hashes.insert(neighbour_hash);
                        node_cache[neighbour_hash].nn_hashes.insert(hash);
                    }
                }
            }
            for (Index idx = 0; auto& nodeEl: node_cache) {
                nodeEl.second.id = idx;
                ++idx;
            }

        }

    }
};

template<floating_point_number Real, integer_number Index>
class PlanarTriangulation{
    Index n_length;
public:
    std::vector<std::vector<Index>> nn_ids;
    std::vector<bool> is_bulk;
    [[nodiscard]] Index ij_to_id(Index i, Index j){return i*n_length+j;}
    [[nodiscard]] Index id_to_i(Index id){return id/n_length;}
    [[nodiscard]] Index id_to_j(Index id){return id%n_length;}

    // TL T TR
    //  L    R
    // BL B BR

    [[nodiscard]] Index TL(Index id){ return ij_to_id(id_to_i(id)-1, id_to_j(id)-1);}
    [[nodiscard]] Index T (Index id){ return ij_to_id(id_to_i(id)-1, id_to_j(id)  );}
    [[nodiscard]] Index TR(Index id){ return ij_to_id(id_to_i(id)-1, id_to_j(id)+1);}
    [[nodiscard]] Index  L(Index id){ return ij_to_id(id_to_i(id)  , id_to_j(id)-1);}
    [[nodiscard]] Index  R(Index id){ return ij_to_id(id_to_i(id)  , id_to_j(id)+1);}
    [[nodiscard]] Index BL(Index id){ return ij_to_id(id_to_i(id)+1, id_to_j(id)-1);}
    [[nodiscard]] Index B (Index id){ return ij_to_id(id_to_i(id)+1, id_to_j(id)  );}
    [[nodiscard]] Index BR(Index id){ return ij_to_id(id_to_i(id)+1, id_to_j(id)+1);}

    [[nodiscard]] std::vector<Index> bulk_odd_j_neighbor_ids(Index id){
        return { B(id), R(id), TR(id), T(id), TL(id), L(id) };
    }

    [[nodiscard]] std::vector<Index> bulk_even_j_neighbor_ids(Index id){
        return { T(id), L(id), BL(id), B(id), BR(id), R(id) };
    }

    [[nodiscard]] std::vector<Index> top_boundary_odd_j_neighbor_ids(Index id){
        return { L(id), B(id), R(id) };
    }

    [[nodiscard]] std::vector<Index> top_boundary_even_j_neighbor_ids(Index id){
        return { L(id), BL(id), B(id), BR(id), R(id) };
    }

    [[nodiscard]] std::vector<Index> bottom_boundary_odd_j_neighbor_ids(Index id){
        return { R(id), TR(id), T(id), TL(id), L(id) };
    }

    [[nodiscard]] std::vector<Index> bottom_boundary_even_j_neighbor_ids(Index id){
        return { T(id), L(id), R(id) };
    }

    [[nodiscard]] std::vector<Index> left_boundary_neighbor_ids(Index id){
        return { T(id), B(id), BR(id), R(id)};
    }

    [[nodiscard]] std::vector<Index> right_boundary_odd_j_neighbor_ids(Index id){
        return { T(id), TL(id), L(id), B(id) };
    }

    [[nodiscard]] std::vector<Index> right_boundary_even_j_neighbor_ids(Index id){
        return { T(id), L(id), BL(id), B(id) };
    }

    PlanarTriangulation(Index n_length_inp, Index n_width):n_length(n_length_inp){
        Index N_nodes = n_length*n_width;

        nn_ids.resize(N_nodes);
        is_bulk.resize(N_nodes,false);
        for(Index i=0; i<n_width;++i) {
            for (Index j = 0; j<n_length; ++j) {
                Index id = ij_to_id(i, j);
            }
        }
        // populate_bulk
        for(Index i=1; i<n_width-1;++i){
            for(Index j=1; j<n_length-1;++j){
                Index bulk_id = ij_to_id(i,j);
                is_bulk[bulk_id] = true;
                if(j%2==0){
                    nn_ids[bulk_id] = bulk_even_j_neighbor_ids(bulk_id);
                }else{
                    nn_ids[bulk_id] = bulk_odd_j_neighbor_ids(bulk_id);
                }
            }
        }

        // populate top and bottom boundaries
        for (Index j = 1; j<n_length-1; ++j) {
            Index i = 0;
            Index id = ij_to_id(i,j);
            if (j%2==0) {
                nn_ids[id] = top_boundary_even_j_neighbor_ids(id);
            }else{
                nn_ids[id] = top_boundary_odd_j_neighbor_ids(id);
            }

            i = n_width-1;
            id = ij_to_id(i,j);
            if (j%2==0) {
                nn_ids[id] = bottom_boundary_even_j_neighbor_ids(id);
            }else{
                nn_ids[id] = bottom_boundary_odd_j_neighbor_ids(id);
            }
        }

        // populate left and right boundaries
        for (Index i = 1; i<n_width-1; ++i) {
            Index j = 0;
            Index id = ij_to_id(i,j);
            nn_ids[id] = left_boundary_neighbor_ids(id);

            j = n_length-1;
            id = ij_to_id(i,j);
            if (j%2==0) {
                nn_ids[id] = right_boundary_even_j_neighbor_ids(id);
            }else{
                nn_ids[id] = right_boundary_odd_j_neighbor_ids(id);
            }
        }

        // populate top left corner
        nn_ids[0] = std::vector<Index>{B(0), BR(0), R(0)};
        // populate bottom left corner
        Index bottom_left_id = ij_to_id(n_width-1,0);
        nn_ids[bottom_left_id] = std::vector<Index>{R(bottom_left_id), T(bottom_left_id)}; //Todo this bond will never flip
        // populate top and bottom right corner
        Index top_right_id = n_length-1;
        Index bottom_right_id = N_nodes-1;
        if((n_length-1)%2==0){
            nn_ids[top_right_id] = {L(top_right_id), BL(top_right_id), B(top_right_id)};
            nn_ids[bottom_right_id] = {T(bottom_right_id), L(bottom_right_id)};// Todo this bond will never flip
        }else{
            nn_ids[top_right_id] = {L(top_right_id), /*BL(top_right_id),*/ B(top_right_id)}; // Todo this bond will never flip
            nn_ids[bottom_right_id] = {T(bottom_right_id), TL(bottom_right_id), L(bottom_right_id),};// Todo this bond will never flip
        }


    }
};
}
#endif //FLIPPY_TRIANGULATOR_HPP
