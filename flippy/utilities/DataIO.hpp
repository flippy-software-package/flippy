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

#ifndef FLIPPY_HPP_DATASAVER_H
#define FLIPPY_HPP_DATASAVER_H

#include "Triangulation.hpp"
#include "custom_concepts.hpp"
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <utility>

namespace fp::experimental {
class DataSaver {
    public:
    virtual void reset_memory()        = 0;
    virtual void stream_data_to_file() = 0;
    void         flush_data_to_file() {
        stream_data_to_file();
        reset_memory();
    }
};

class jsonDataSaver : public DataSaver {
    std::string   data_;
    std::ofstream data_file_;

    public:
    explicit jsonDataSaver(const std::filesystem::path &save_file_path) : data_(), data_file_(save_file_path) {}

    void reset_memory() override { data_.clear(); }
    void save_current_state_in_memory(const Triangulation &triangulation) {
        data_ += triangulation.nodes().make_data().dump();
    }
    void stream_data_to_file() override { data_file_ << data_; }
};

enum struct XYZ_TYPE : std::uint8_t { S = 0, R, I };

static std::string xyz_to_string(XYZ_TYPE const &xyz_type) noexcept(false) {
    switch (xyz_type) {
        using enum XYZ_TYPE;
    case S:
        return "S";
    case R:
        return "R";
    case I:
        return "I";
    default:
        throw std::runtime_error("XYZ_TYPE not recognized, it must be S, R or I");
    }
}

struct xyzProperty {
    std::string name;
    XYZ_TYPE    xyz_type;
    short       column_count;
};

class xyzDataSaver : public DataSaver {
    std::string   data_;
    std::ofstream data_file_;

    std::string size_string{};
    std::string properties_string;

    using StreamParticleFuncType = std::function<std::string(const Node &, const Triangulation &)>;
    StreamParticleFuncType stream_particle;
    using GlobalsDictType      = std::unordered_map<std::string, std::string>;
    using GlobalsMakerFuncType = std::function<GlobalsDictType(const Triangulation &)>;

    void set_up_property_string(const std::vector<xyzProperty> &properties) {
        properties_string.append("Lattice=\"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\" Properties=");
        for (const auto &property : properties) {
            properties_string.append(property.name + ":" + xyz_to_string(property.xyz_type) + ":" +
                                     std::to_string(property.column_count) + ":");
        }
        properties_string.pop_back();
    }

    void make_header(const Triangulation &triangulation) {
        size_string = std::to_string(triangulation.size());
        data_.append(size_string);
        data_.append("\n");
        data_.append(properties_string);
    }

    void make_globals(const GlobalsDictType &globals_dict) {
        std::string prop{};
        prop.reserve(200);
        for (auto &[key, val] : globals_dict) {
            prop.append(" ");
            prop.append(key);
            prop.append("=");
            prop.append(val);
        }
        data_.append(prop);
    }
    void save_nodes(const Triangulation &triangulation) {
        data_.append("\n");
        for (const auto &node : triangulation.nodes()) {
            data_.append(stream_particle(node, triangulation));
        }
    }

    public:
    explicit xyzDataSaver(const std::filesystem::path    &save_file_path,
                          const std::vector<xyzProperty> &properties_inp,
                          StreamParticleFuncType          stream_particle_inp)
        : data_(), data_file_(save_file_path), stream_particle(std::move(stream_particle_inp)) {
        set_up_property_string(properties_inp);
    }

    void reset_memory() override { data_.clear(); }

    void save_current_state_in_memory(const Triangulation &triangulation) {
        make_header(triangulation);
        save_nodes(triangulation);
    }

    void save_current_state_in_memory(const Triangulation &triangulation, const GlobalsDictType &globals_dict) {
        make_header(triangulation);
        make_globals(globals_dict);
        save_nodes(triangulation);
    }

    void save_current_state_in_memory(const Triangulation &triangulation, const GlobalsMakerFuncType &globals_maker) {
        make_header(triangulation);
        auto globals_dict = globals_maker(triangulation);
        make_globals(globals_dict);
        save_nodes(triangulation);
    }

    void stream_data_to_file() override { data_file_ << data_; }
};

class csvGlobalsSaver : public DataSaver {
    std::ofstream            file;
    std::string              data;
    std::string              filename;
    std::vector<std::string> names;
    std::vector<std::string> values;
    std::string              header;
    std::string              sep;
    bool                     header_written = false;

    public:
    csvGlobalsSaver(std::string filename_inp, std::vector<std::string> names_inp, std::string sep_inp = ",")
        : filename(std::move(filename_inp)), names(std::move(names_inp)), sep(std::move(sep_inp)) {
        file.open(filename);
        header = names[0];
        for (size_t i = 1; i < names.size(); ++i) {
            header += sep + names[i];
        }
        header_written = false;
    }

    void reset_memory() override { data.clear(); }

    void save_current_state_in_memory(const std::vector<double> &values_inp) {
        if (!header_written) {
            data.append(header).append("\n");
            header_written = true;
        }
        std::string row = std::to_string(values_inp[0]);
        for (size_t i = 1; i < values_inp.size(); ++i) {
            row += sep + std::to_string(values_inp[i]);
        }
        data.append(row).append("\n");
    }

    void stream_data_to_file() override { file << data; }

    ~csvGlobalsSaver() { file.close(); }
};

} // namespace fp::experimental
#endif
