#ifndef FLIPPY_TRIANGULATOR_HPP
#define FLIPPY_TRIANGULATOR_HPP

#include <array>
#include <vector>
#include <unordered_set>
#include "vec3.hpp"

namespace fp::implementation{

template<typename Real, typename Index>
struct SimpleNodeData{
  std::string hash{};
  Index id{};
  vec3<Real> pos{};
  std::unordered_set<std::string> nn_hashes{};
};

template<typename Index> requires std::is_integral_v<Index> std::string hash_node(Index i){ return std::to_string(i);}

template<typename Index> requires std::is_integral_v<Index> std::string hash_node(Index i, Index j, Index n){
    Index a = std::min(i, j);
    Index b = std::max(i, j);
    return std::to_string(a)+"_"+ std::to_string(b) + "_" + std::to_string(n);
}

template<typename Index> requires std::is_integral_v<Index> std::string hash_node(Index c0, Index c1, Index c2,
        Index i, Index j){
    std::vector<Index> cv{c0, c1, c2};
    std::sort(cv.begin(), cv.end());

    return std::to_string(cv[0]) + "_" + std::to_string(cv[1]) + "_" + std::to_string(cv[2])
    + "_" + std::to_string(i) + "_" + std::to_string(j);
}

constexpr int N_ICOSA_FACEs=20;
constexpr int N_ICOSA_EDGEs=30;
constexpr int N_ICOSA_NODEs=12;

//constexpr std::array<std::array<int,2>,N_ICOSA_EDGEs> ICOSAHEDRON_EDGES = {
//        std::array<int,2>{0, 2},
//        std::array<int,2>{0, 4},
//        std::array<int,2>{0, 5},
//        std::array<int,2>{0, 8},
//        std::array<int,2>{0, 9},
//        std::array<int,2>{1, 3},
//        std::array<int,2>{1, 6},
//        std::array<int,2>{1, 7},
//        std::array<int,2>{1,10},
//        std::array<int,2>{1,11},
//        std::array<int,2>{2, 6},
//        std::array<int,2>{2, 7},
//        std::array<int,2>{2, 8},
//        std::array<int,2>{2, 9},
//        std::array<int,2>{3, 4},
//        std::array<int,2>{3, 5},
//        std::array<int,2>{3,10},
//        std::array<int,2>{3, 11},
//        std::array<int,2>{4, 5},
//        std::array<int,2>{4, 8},
//        std::array<int,2>{4, 10},
//        std::array<int,2>{5, 9},
//        std::array<int,2>{5, 11},
//        std::array<int,2>{6, 7},
//        std::array<int,2>{6, 8},
//        std::array<int,2>{6, 10},
//        std::array<int,2>{7, 9},
//        std::array<int,2>{7, 11},
//        std::array<int,2>{8, 10},
//        std::array<int,2>{9, 11}
//};

constexpr std::array<int, N_ICOSA_FACEs> FACE_IDs{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
//constexpr std::array<int, N_ICOSA_FACEs> FACE_IDs{0, 1, 2, 3, 4, 5, 14, 13, 12, 11, 10, 9, 8, 7, 6, 16, 15, 19, 18, 17};
constexpr std::array<int, N_ICOSA_NODEs> NODE_IDs{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
constexpr std::array<std::array<int,3>,N_ICOSA_FACEs> FACE_CORNER_NODES_={
        std::array<int,3>{NODE_IDs[0 ], NODE_IDs[5 ], NODE_IDs[1 ]},
        std::array<int,3>{NODE_IDs[0 ], NODE_IDs[1 ], NODE_IDs[2 ]},
        std::array<int,3>{NODE_IDs[0 ], NODE_IDs[2 ], NODE_IDs[3 ]},
        std::array<int,3>{NODE_IDs[0 ], NODE_IDs[3 ], NODE_IDs[4 ]},
        std::array<int,3>{NODE_IDs[0 ], NODE_IDs[4 ], NODE_IDs[5 ]},
        std::array<int,3>{NODE_IDs[10], NODE_IDs[5 ], NODE_IDs[4 ]},
        std::array<int,3>{NODE_IDs[10], NODE_IDs[6 ], NODE_IDs[5 ]},
        std::array<int,3>{NODE_IDs[6 ], NODE_IDs[1 ], NODE_IDs[5 ]},
        std::array<int,3>{NODE_IDs[6 ], NODE_IDs[7 ], NODE_IDs[1 ]},
        std::array<int,3>{NODE_IDs[7 ], NODE_IDs[2 ], NODE_IDs[1 ]},
        std::array<int,3>{NODE_IDs[7 ], NODE_IDs[8 ], NODE_IDs[2 ]},
        std::array<int,3>{NODE_IDs[8 ], NODE_IDs[3 ], NODE_IDs[2 ]},
        std::array<int,3>{NODE_IDs[8 ], NODE_IDs[9 ], NODE_IDs[3 ]},
        std::array<int,3>{NODE_IDs[9 ], NODE_IDs[4 ], NODE_IDs[3 ]},
        std::array<int,3>{NODE_IDs[9 ], NODE_IDs[10], NODE_IDs[4 ]},
        std::array<int,3>{NODE_IDs[11], NODE_IDs[10], NODE_IDs[9 ]},
        std::array<int,3>{NODE_IDs[11], NODE_IDs[6 ], NODE_IDs[10]},
        std::array<int,3>{NODE_IDs[11], NODE_IDs[7 ], NODE_IDs[6 ]},
        std::array<int,3>{NODE_IDs[11], NODE_IDs[8 ], NODE_IDs[7 ]},
        std::array<int,3>{NODE_IDs[11], NODE_IDs[9 ], NODE_IDs[8 ]}
};

//const std::array<std::array<int,3>,N_ICOSA_FACEs> FACE_CORNER_NODES={
//        FACE_CORNER_NODES_[FACE_IDs[0]],
//        FACE_CORNER_NODES_[FACE_IDs[1]],
//        FACE_CORNER_NODES_[FACE_IDs[2]],
//        FACE_CORNER_NODES_[FACE_IDs[3]],
//        FACE_CORNER_NODES_[FACE_IDs[4]],
//        FACE_CORNER_NODES_[FACE_IDs[5]],
//        FACE_CORNER_NODES_[FACE_IDs[6]],
//        FACE_CORNER_NODES_[FACE_IDs[7]],
//        FACE_CORNER_NODES_[FACE_IDs[8]],
//        FACE_CORNER_NODES_[FACE_IDs[9]],
//        FACE_CORNER_NODES_[FACE_IDs[10]],
//        FACE_CORNER_NODES_[FACE_IDs[11]],
//        FACE_CORNER_NODES_[FACE_IDs[12]],
//        FACE_CORNER_NODES_[FACE_IDs[13]],
//        FACE_CORNER_NODES_[FACE_IDs[14]],
//        FACE_CORNER_NODES_[FACE_IDs[15]],
//        FACE_CORNER_NODES_[FACE_IDs[16]],
//        FACE_CORNER_NODES_[FACE_IDs[17]],
//        FACE_CORNER_NODES_[FACE_IDs[18]],
//        FACE_CORNER_NODES_[FACE_IDs[19]]
//};
const std::array<std::array<int,3>,N_ICOSA_FACEs> FACE_CORNER_NODES={
        std::array<int,3>{0, 5, 1},
        std::array<int,3>{0, 1, 2},
        std::array<int,3>{0, 2, 3},
        std::array<int,3>{0, 3, 4},
        std::array<int,3>{0, 4, 5},
        std::array<int,3>{10, 5, 4},
        std::array<int,3>{10, 6, 5},
        std::array<int,3>{6, 1, 5},
        std::array<int,3>{6, 7, 1},
        std::array<int,3>{7, 2, 1},
        std::array<int,3>{7, 8, 2},
        std::array<int,3>{8, 3, 2},
        std::array<int,3>{8, 9, 3},
        std::array<int,3>{9, 4, 3},
        std::array<int,3>{9, 10, 4},
        std::array<int,3>{11, 10, 9},
        std::array<int,3>{11, 6, 10},
        std::array<int,3>{11, 7, 6},
        std::array<int,3>{11, 8, 7},
        std::array<int,3>{11, 9, 8}
};


//template<typename Index>
//requires std::is_integral_v<Index>
//Index k_lt(Index i, Index j){return (i*(i+1)>>1) +j;}
//
//template<typename Index>
//requires std::is_integral_v<Index>
//Index i_lt(Index k){return floor(0.5f*(sqrt(1+8*k)-1.f));}
//
//template<typename Index>
//requires std::is_integral_v<Index>
//Index j_lt(Index k, Index i){return k - (i*(i+1)>>1);}
//
//template<typename Index>
//requires std::is_integral_v<Index>
//std::pair<Index, Index> ij_lt(Index k){
//    Index i = i_lt(k);
//    return {i, j_lt(k, i)};
//}

template<typename Real>
vec3<Real> r_S1(Real R, Real t, Real f) {
    vec3<Real> r{R * std::sin(t) * std::cos(f), R * std::sin(t) * std::sin(f), R * std::cos(t)};
    return r;
}

template<typename Real, typename Index>
std::unordered_map<std::string, SimpleNodeData<Real, Index>> make_corner_nodes() {

    Real R = 1.;
    std::unordered_map<std::string, SimpleNodeData<Real, Index>> base_nodes(N_ICOSA_NODEs);
    base_nodes[hash_node(0)] = {.hash=hash_node(0), .id=0, .pos=r_S1<Real>(R, 0., 0.)};
    std::string hash;
    hash.reserve(2);
    for (Index i = 1; i < 6; ++i) {
        hash = hash_node(i);
        base_nodes[hash] = {
                .hash=hash,
                .id=i,
                .pos=r_S1<Real>(R, M_PI/2. - std::atan(0.5), 2*M_PI*(i-1.)/5.)};
    }

    for (Index i = 6; i < N_ICOSA_NODEs - 1; ++i) {
        hash = hash_node(i);
        base_nodes[hash] = {
                .hash=hash,
                .id=i,
                .pos=r_S1<Real>(R, M_PI/2. + std::atan(0.5), 2*M_PI*(i-6.5)/5.)};
    }
    hash = hash_node(N_ICOSA_NODEs-1);
    base_nodes[hash] = {
            .hash=hash,
            .id=static_cast<Index>(N_ICOSA_NODEs-1),
            .pos=r_S1<Real>(R, M_PI, 0.)};
    return base_nodes;
}

enum TriangleRegion {TOP_CORNER, BOTTOM_LEFT_CORNER, BOTTOM_RIGHT_CORNER, LEFT_EDGE, BOTTOM_EDGE, DIAGONAL_EDGE, BULK};

template<typename Index>
TriangleRegion get_region(Index i, Index j, Index sizeMinOne){
    if(i==0){return TOP_CORNER;}
    else if(j==0 && i==sizeMinOne){return BOTTOM_LEFT_CORNER;}
    else if(j==sizeMinOne && i==sizeMinOne){return BOTTOM_RIGHT_CORNER;}
    else if(j==0){return LEFT_EDGE;}
    else if(i==sizeMinOne){return BOTTOM_EDGE;}
    else if(i==j){return DIAGONAL_EDGE;}
    else{return BULK;}
}

//std::string hash_any(Index c0, Index c1, Index c2, Index i, Index j, Index maxIdx){
template<typename Index>
std::string hash_any(Index c0, Index c1, Index c2, Index i, Index j, Index maxIdx){
    switch (get_region(i, j, maxIdx)) {
        case TOP_CORNER:
            return hash_node(c0);
        case BOTTOM_LEFT_CORNER:
            return hash_node(c1);
        case BOTTOM_RIGHT_CORNER:
            return hash_node(c2);
        case LEFT_EDGE:
            return hash_node(c0,c1,i);
        case BOTTOM_EDGE:
            return hash_node(c1,c2,j);
        case DIAGONAL_EDGE:
            return hash_node(c0,c2,j);
        case BULK:
            return hash_node(c0,c1,c2, i, j);
    }
}

template<typename Index>
std::vector<std::string> neighbour_hash_vec(Index c0, Index c1, Index c2, Index i, Index j, Index maxIdx){
    std::vector<std::string>  neighbour_hash; neighbour_hash.reserve(6);
    switch (get_region(i, j, maxIdx)){

    case TOP_CORNER:
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,1,0,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,1,1,maxIdx));
        return neighbour_hash;

    case BOTTOM_LEFT_CORNER:
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i,j+1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i-1,j,maxIdx));
        return neighbour_hash;

    case BOTTOM_RIGHT_CORNER:
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i,j-1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i-1,j-1,maxIdx));
        return neighbour_hash;

    case LEFT_EDGE:
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i-1, j,   maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i+1, j,   maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i,   j+1, maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i+1, j+1, maxIdx));
        return neighbour_hash;

    case BOTTOM_EDGE:
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i,   j-1,   maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i-1, j-1, maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i-1, j,   maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i,  j+1,  maxIdx));
        return neighbour_hash;

    case DIAGONAL_EDGE:
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i, j-1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i-1, j-1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i+1, j,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i+1, j+1,maxIdx));
        return neighbour_hash;

    case BULK:
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i, j-1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i-1, j-1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2,i-1, j,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i, j+1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i+1, j+1,maxIdx));
        neighbour_hash.push_back(hash_any<Index>(c0,c1,c2, i+1, j,maxIdx));
        return neighbour_hash;
}

}

template<typename Real, typename Index>
vec3<Real> get_pos(vec3<Real>const& p0, vec3<Real>const& p1,
                   vec3<Real>const& p2, Index i, Index j, Index maxIdx){

    Real wi = static_cast<Real>(i)/(static_cast<Real>(maxIdx));
    Real wj = static_cast<Real>(j)/(static_cast<Real>(maxIdx));
    switch (get_region(i, j, maxIdx)) {
        case TOP_CORNER:
            return p0;
        case BOTTOM_LEFT_CORNER:
            return p1;
        case BOTTOM_RIGHT_CORNER:
            return p2;
        case LEFT_EDGE:
             return p0 + wi*(p1 - p0);
        case BOTTOM_EDGE:
          return p1 + wj*(p2 - p1);
        case DIAGONAL_EDGE:
             return p0 + wj*(p2 - p0);
        case BULK:
            return p0 + wi*(p1 - p0) + wj*(p2 - p1);
//            return p0 + wi*(p1 - p0)/(p1 - p0).norm() + wj*(p2 - p0)/(p2 - p0).norm();
    }
}

template<typename Index>
std::tuple<Index, Index, Index> get_sorted_face_nodes(std::array<int,3> face){
    std::sort(face.begin(), face.end());
    return {static_cast<Index>(face[0]),
            static_cast<Index>(face[1]),
            static_cast<Index>(face[2])};
}

template<typename Real, typename Index>
void make_face_nodes(std::unordered_map<std::string,SimpleNodeData<Real, Index>>& node_cache, Index nIter){
//    vec3<Real> e1, e2, anchor, pos;
    vec3<Real> p0, p1, p2, pos;
    Index nEdge = nIter +2;// total Number of nodes on an edge
    Index maxIdx = nIter +1;// max value i or j can have

    std::string hash;
    hash.reserve(10);

    std::string c0_h, c1_h, c2_h;
    for (auto face : FACE_CORNER_NODES) {
        auto [c0, c1, c2] = get_sorted_face_nodes<Index>(face);
        c0_h = hash_node(c0); c1_h = hash_node(c1); c2_h = hash_node(c2);
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

    for (auto face : FACE_CORNER_NODES) {
//        std::sort(face.begin(),face.end());
//        auto [c0, c1, c2] = face;
        auto [c0, c1, c2] = get_sorted_face_nodes<Index>(face);
        for (Index i = 0; i<nEdge; ++i) {
            for (Index j = 0; j<=i; ++j) {
                hash = hash_any(c0, c1, c2, i, j, maxIdx);
                std::vector<std::string> neighbour_hashes = neighbour_hash_vec(c0, c1, c2, i, j, maxIdx);
                for(auto const& neighbour_hash: neighbour_hashes){
                    node_cache[hash].nn_hashes.insert(neighbour_hash);
                    node_cache[neighbour_hash].nn_hashes.insert(hash);
                }
            }
        }
        for(Index idx=0; auto & nodeEl :node_cache){
            nodeEl.second.id = idx;
            ++idx;
        }

    }

}
}
#endif //FLIPPY_TRIANGULATOR_HPP
