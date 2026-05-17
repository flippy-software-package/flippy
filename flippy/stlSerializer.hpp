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

#ifndef FLIPPY_HPP_STLSERIALIZER_HPP
#define FLIPPY_HPP_STLSERIALIZER_HPP

#include "vec3.hpp"
#include <cstdint>
#include <fstream>
#include <optional>
#include <vector>

namespace fp::implementation {
struct rawSTLTriangle {
    float normal[3]{};
    float v1[3]{};
    float v2[3]{};
    float v3[3]{};

    uint16_t attr_byte_count{0};
};

struct rawSTLSolid {
    char                        header[80]{};
    std::uint32_t               num_triangles{};
    std::vector<rawSTLTriangle> stl_triangles;
};

[[maybe_unused]] static void printSTLTriangle(const rawSTLTriangle             &triangle,
                                              const std::optional<std::string> &name = std::nullopt) {
    if (name) { std::cout << name.value() << '\n'; }
    std::cout << triangle.normal[0] << ' ' << triangle.normal[1] << ' ' << triangle.normal[2] << '\n';
    std::cout << triangle.v1[0] << ' ' << triangle.v1[1] << ' ' << triangle.v1[2] << '\n';
    std::cout << triangle.v2[0] << ' ' << triangle.v2[1] << ' ' << triangle.v2[2] << '\n';
    std::cout << triangle.v3[0] << ' ' << triangle.v3[1] << ' ' << triangle.v3[2] << '\n';
    std::cout << "------------------------------------------------------------------\n";
}

static inline rawSTLTriangle loadSTLTriangle(std::ifstream &inp) {
    rawSTLTriangle triangle;
    inp.read(reinterpret_cast<char *>(&triangle), 50 /*sizeof(STLTriangle)*/);
    return triangle;
}

static inline rawSTLSolid loadSTLSolid(std::ifstream &inp) {
    rawSTLSolid solid;
    inp.read(reinterpret_cast<char *>(&solid.header), sizeof(rawSTLSolid::header));
    inp.read(reinterpret_cast<char *>(&solid.num_triangles), sizeof(rawSTLSolid::num_triangles));
    std::vector<rawSTLTriangle> triangles;
    triangles.reserve(solid.num_triangles);
    //        std::cout << "loading STL solid with: " << solid.num_triangles << " triangles'\n";
    for (uint32_t i = 0; i < solid.num_triangles; ++i) {
        triangles.push_back(loadSTLTriangle(inp));
    }
    solid.stl_triangles = std::move(triangles);
    return solid;
}

template <fp::floating_point_number Real, fp::IndexingNumber Index> struct stlNode {
    Index      id;
    vec3<Real> pos;

    bool operator==(const stlNode &other) const { return pos == other.pos; }
};

template <fp::floating_point_number Real, fp::IndexingNumber Index> class stlTriangle {
    vec3<Real> normal_{};

    static vec3<Real>
    make_normal(const stlNode<Real, Index> &n1, const stlNode<Real, Index> &n2, const stlNode<Real, Index> &n3) {
        return -1 * (n2.pos - n1.pos).cross(n3.pos - n1.pos).normalize();
    }

    public:
    stlNode<Real, Index> n1_;
    stlNode<Real, Index> n2_;
    stlNode<Real, Index> n3_;

    stlTriangle(stlNode<Real, Index> &n1_inp,
                stlNode<Real, Index> &n2_inp,
                stlNode<Real, Index> &n3_inp,
                const vec3<float>    &normal_inp)
        : normal_(normal_inp), n1_(n1_inp), n2_(n2_inp), n3_(n3_inp) {}

    stlTriangle(stlNode<Real, Index> &n1_inp, stlNode<Real, Index> &n2_inp, stlNode<Real, Index> &n3_inp)
        : stlTriangle(n1_inp, n2_inp, n3_inp, make_normal(n1_inp, n2_inp, n3_inp)) {}

    explicit stlTriangle(const rawSTLTriangle &stl_triangle)
        : normal_({stl_triangle.normal[0], stl_triangle.normal[1], stl_triangle.normal[2]}),
          n1_({.id = 0, .pos = {stl_triangle.v1[0], stl_triangle.v1[1], stl_triangle.v1[2]}}),
          n2_({.id = 0, .pos = {stl_triangle.v2[0], stl_triangle.v2[1], stl_triangle.v2[2]}}),
          n3_({.id = 0, .pos = {stl_triangle.v3[0], stl_triangle.v3[1], stl_triangle.v3[2]}}) {}

    [[nodiscard]] const vec3<float> &normal() const { return normal_; }

    [[nodiscard]] rawSTLTriangle to_stl_triangle() const {

        rawSTLTriangle tr{.normal          = {normal_.x, normal_.y, normal_.z},
                          .v1              = {n1_.pos.x, n1_.pos.y, n1_.pos.z},
                          .v2              = {n2_.pos.x, n2_.pos.y, n2_.pos.z},
                          .v3              = {n3_.pos.x, n3_.pos.y, n3_.pos.z},
                          .attr_byte_count = 0};
        return tr;
    };

    std::pair<Index, Index> other_node_ids(Index loc_node_id) {
        switch (loc_node_id) {
        case 0:
            return {n2_.id, n3_.id};
        case 1:
            return {n1_.id, n3_.id};
        case 2:
            return {n1_.id, n2_.id};
        default:
            std::cerr << loc_node_id << "is out of range for as stlTriangl index";
            exit(12);
        }
    }

    stlNode<Real, Index> &operator[](Index idx) {
        switch (idx) {
        case 0:
            return n1_;
        case 1:
            return n2_;
        case 2:
            return n3_;
        default:
            std::cerr << idx << "is out of range for as stlTriangl index";
            exit(12);
        }
    }
};

template <fp::floating_point_number Real, fp::IndexingNumber Index> class stlSerializer {
    static std::vector<stlTriangle<Real, Index>> solidTostlTriangls(const rawSTLSolid &solid) {
        std::vector<stlTriangle<Real, Index>> triangles;
        triangles.reserve(solid.num_triangles);
        for (const auto &stl_triangle : solid.stl_triangles) {
            triangles.emplace_back(stl_triangle);
        }
        return triangles;
    }

    public:
    static std::vector<stlTriangle<Real, Index>>
    read_STLSolid_into_triangle_vec(const std::filesystem::path &input_file) {
        std::ifstream inp(input_file, std::ios::in | std::ios::binary);
        rawSTLSolid   solid{loadSTLSolid(inp)};
        inp.close();

        std::vector<stlTriangle<Real, Index>> triangles{solidTostlTriangls(solid)};

        return triangles;
    }
};

} // namespace fp::implementation

#endif // FLIPPY_HPP_STLSERIALIZER_HPP
