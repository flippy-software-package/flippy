#ifndef FLIPPY_HPP_DATASAVER_H
#define FLIPPY_HPP_DATASAVER_H

#include "custom_concepts.hpp"
#include "Triangulation.hpp"
#include <filesystem>
#include <iostream>
#include <utility>

namespace fp::experimental{
    class DataSaver {
    public:
        virtual void reset_memory() = 0;
        virtual void stream_data_to_file() = 0;
        void flush_data_to_file(){
            stream_data_to_file();
            reset_memory();
        }

    };

    class jsonDataSaver: public DataSaver{
        std::string data_;
        std::ofstream data_file_;
    public:

        explicit jsonDataSaver(std::filesystem::path const& save_file_path): data_(), data_file_(save_file_path) { }

        void reset_memory() override { data_.clear(); }
        void save_current_state_in_memory(Triangulation const & triangulation) {

            data_ += triangulation.nodes().make_data().dump();
        }
        void stream_data_to_file() override { data_file_ << data_; }
    };

    enum XYZ_TYPE {S, R, I};

    static std::string xyz_to_string(XYZ_TYPE const& xyz_type) noexcept(false){
        switch(xyz_type){
            case S: return "S";
            case R: return "R";
            case I: return "I";
            default: throw std::runtime_error("XYZ_TYPE not recognized, it must be S, R or I");
        }
    }

    struct xyzProperty
    {
        std::string name;
        XYZ_TYPE xyz_type;
        short column_count;
    };

    class xyzDataSaver :public DataSaver{
        std::string data_;
        std::ofstream data_file_;

        std::string size_string{};
        std::string properties_string;

        using StreamParticleFuncType = std::function<std::string(Node const&, Triangulation const&)>;
        StreamParticleFuncType stream_particle;
        using GlobalsMakerFuncType = std::function<std::unordered_map<std::string, std::string>(Triangulation const&)>;

        void set_up_property_string(std::vector<xyzProperty> const& properties)
        {
            properties_string.append("Lattice=\"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\" Properties=");
            for (auto const &property : properties) {
                properties_string.append(property.name + ":" + xyz_to_string(property.xyz_type) + ":" + std::to_string(property.column_count) + ":");
            }
            properties_string.pop_back();
            properties_string.append("\n");
        }
    public:

        explicit xyzDataSaver(std::filesystem::path const& save_file_path,
                              std::vector<xyzProperty> const& properties_inp,
                              StreamParticleFuncType  stream_particle_inp
        ): data_(), data_file_(save_file_path), stream_particle(std::move(stream_particle_inp)){
            set_up_property_string(properties_inp);
        }

        void reset_memory() override { data_.clear(); }

        void save_current_state_in_memory(Triangulation const & triangulation,
                                          GlobalsMakerFuncType const& globals_maker ) {
            size_string = std::to_string(triangulation.size());
            data_.append(size_string);
            data_.append("\n");
            std::string prop{properties_string};
            prop.reserve(500);
            for(auto & [key, val]: globals_maker(triangulation)){
                prop.append(" ");
                prop.append(key);
                prop.append("=");
                prop.append(val);
            }
            data_.append(prop);
            data_.append("\n");

            for(const auto& node: triangulation.nodes()){ data_.append(stream_particle(node, triangulation)); }
        }

        void stream_data_to_file() override { data_file_ << data_; }
    };

    class csvGlobalsSaver: public DataSaver
    {
        std::ofstream file;
        std::string data;
        std::string filename;
        std::vector<std::string> names;
        std::vector<std::string> values;
        std::string header;
        std::string sep;
        bool header_written = false;

    public:
        csvGlobalsSaver(std::string filename_inp, std::vector<std::string> names_inp, std::string sep_inp=",")
                : filename(std::move(filename_inp)), names(std::move(names_inp)), sep(std::move(sep_inp))
        {
            file.open(filename);
            header = names[0];
            for(size_t i=1; i<names.size(); ++i){
                header += sep + names[i];
            }
            header_written = false;
        }

        void reset_memory() override { data.clear(); }

        void save_current_state_in_memory(std::vector<double>const& values_inp){
            if(!header_written){
                data.append(header).append("\n");
                header_written = true;
            }
            std::string row = std::to_string(values_inp[0]);
            for(size_t i=1; i<values_inp.size(); ++i){
                row += sep + std::to_string(values_inp[i]);
            }
            data.append(row).append("\n");
        }

        void stream_data_to_file() override { file << data; }

        ~csvGlobalsSaver(){
            file.close();
        }
    };

}
#endif
