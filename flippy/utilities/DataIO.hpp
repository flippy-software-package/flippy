#ifndef FLIPPY_HPP_DATASAVER_H
#define FLIPPY_HPP_DATASAVER_H

#include "custom_concepts.hpp"
#include "Triangulation.hpp"
#include <filesystem>
#include <iostream>

namespace fp{
    template <typename DataIO>
    concept DataIOInterface = requires(DataIO D, DataIO const constD, Triangulation<SPHERICAL_TRIANGULATION> const&
            triangulation) {
        typename DataIO::result_type;
        { D.save_current_state_in_memory(triangulation) } -> std::same_as<DataIO>;
        { D.reset_memory() } -> std::same_as<void>;
        { D.stream_data_to_file()} -> std::same_as<void>;
        { D.load_from_file() } -> std::same_as<Triangulation<SPHERICAL_TRIANGULATION>>;
        { D.get_file_handle() } -> std::same_as<std::ofstream const &>;
        { constD.get_file_handle() } -> std::same_as<std::ofstream const &>;
    };

    class xyzDataIO {
        std::string data_;
        std::ofstream data_file_;
    public:

        explicit xyzDataIO(std::filesystem::path const& save_file_path): data_(), data_file_(save_file_path) { }


        void reset_memory() { data_.clear(); }
        template<fp::floating_point_number Real, fp::indexing_number Index>
        void save_current_state_in_memory(Triangulation<SPHERICAL_TRIANGULATION> const & triangulation) { data_ += triangulation.make_egg_data().dump(); }
        void stream_data_to_file(){ data_file_ << data_; }
    };

//    static_assert(DataIOInterface<xyzDataIO, float, unsigned int>, "xyzDataIO is not implementing DataIOInterface");
}
#endif
