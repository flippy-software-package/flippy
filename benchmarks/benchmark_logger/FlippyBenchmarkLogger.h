#ifndef NODES_HPP_FLIPPYBENCHMARKLOGGER_H
#define NODES_HPP_FLIPPYBENCHMARKLOGGER_H
#include <concepts>
#include <flippy.hpp>
#include <code_utils.hpp>

template <typename U>
concept Updater = requires(U u,
                           fp::Triangulation<fp::SPHERICAL_TRIANGULATION> triangulation,
                           fp::Node node, fp::Index nn_id,
                           fp::vec3<fp::Real> displacement,
                           unsigned
                           long long z,
                           std::ostream& os) {
//  { U(triangulation, node) } -> std::same_as<U>;

  {u.move_MC_updater(node, displacement)} -> std::same_as<std::optional<fp::Real>>;
  {u.flip_MC_updater(node, nn_id)} -> std::same_as<void>;

  {u.move_attempt_count()} -> std::same_as<unsigned long>;
  {u.bond_length_move_rejection_count()} -> std::same_as<unsigned long>;
  {u.move_back_count()} -> std::same_as<unsigned long>;

  {u.flip_attempt_count()} -> std::same_as<unsigned long>;
  {u.bond_length_flip_rejection_count()} -> std::same_as<unsigned long>;
  {u.flip_back_count()} -> std::same_as<unsigned long>;

};

template <Updater Updater>
class FlippyBenchmarkLogger {
  Updater const& mcu_;
  fp::Triangulation<fp::SPHERICAL_TRIANGULATION> const& triangulation_;
  cutils::Timer timer{};

FlippyBenchmarkLogger(Updater const& mcu,
                      fp::Triangulation<fp::SPHERICAL_TRIANGULATION> const& triangulation
                        ): mcu_(mcu), triangulation_(triangulation) {};

public:
  static FlippyBenchmarkLogger start_benchmark(Updater const& mcu, fp::Triangulation<fp::SPHERICAL_TRIANGULATION> const& triangulation) {
    FlippyBenchmarkLogger logger(mcu, triangulation);

    PRINT("each of ", triangulation.size(), " beads will be attempted to be moved: ",
          "XXX", "(predicted) times");
    PRINT("leading to a total of", "XXX", "(predicted) MC moves");

    return logger;
  }

std::string log(std::string const& log_dir, fp::Json const& config_data){
  cutils::HumanReadableTime elapsed_time = timer.stop();
  // MonteCarloUpdater counts the number of accepted and rejected moves, distinguishing whether a rejection occurred because of the energy or the bond length constraint.
  // We can use this to print simple statistics here. For example, this will help us decide if our displacement size is too large.
  auto attempt = mcu_.move_attempt_count();
  auto put_back = mcu_.move_back_count();
  auto bond_length_legal_move = attempt - mcu_.bond_length_move_rejection_count();

  auto attempt_flip = mcu_.flip_attempt_count();
  auto flip_back = mcu_.flip_back_count();

  auto topologically_successful_flip = mcu_.flip_attempt_count() - flip_back;
  PRINT("percentage of failed moves: ",
        static_cast<long double>(mcu_.move_back_count() + mcu_.bond_length_move_rejection_count()) / static_cast<long double>(mcu_.move_attempt_count() ));
  PRINT("percentage of failed flips: ",
        static_cast<long double>(mcu_ .flip_back_count() + mcu_.bond_length_flip_rejection_count()) / static_cast<long double>(mcu_.flip_attempt_count()));

  long double time_in_seconds = ((long double)elapsed_time.diff_ns)*1.e-9l;
  long double moves_ps = static_cast<long double>(attempt)/static_cast<long double>(time_in_seconds);
  long double e_moves_ps = static_cast<long double>(bond_length_legal_move)/static_cast<long double>(time_in_seconds);
  long double flips_ps = static_cast<long double>(attempt_flip)/static_cast<long double>(time_in_seconds);
  long double e_flips_ps = static_cast<long double>(topologically_successful_flip)/static_cast<long double>(time_in_seconds);

  auto now = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::tm* time = std::localtime(&now_time);
  std::stringstream date, time_stamp;
  date<<1900+time->tm_year<<"-"<<time->tm_mon+1<<"-"<<time->tm_mday;
  time_stamp<<time->tm_hour<<"-"<<time->tm_min<<"-"<<time->tm_sec;

  std::string log_sub_dir = log_dir + date.str() + "/";
  std::filesystem::create_directories(log_sub_dir);

  std::string log_name = log_sub_dir + time_stamp.str() + ".yml";

  std::ofstream log_file(log_name);
  log_file<<"BENCHMARK VERSION: "<<config_data["benchmark_version"].get<int>()<< '\n';
  log_file<<"config: "<<config_data<<'\n';
  log_file<<"all_attempted_bead_moves: "<<attempt<<'\n';
  log_file<<"expensive_moves: "<< bond_length_legal_move<<'\n';
  log_file<<"rejected_moves_stochastic: "<<put_back<<'\n';
  log_file<<"all_attempted_flips: "<<attempt_flip<<'\n';
  log_file<<"expensive_flips: "<<topologically_successful_flip<<'\n';
  log_file<<"rejected_flip_stochastic: "<<flip_back<<'\n';
  log_file<<"total_simulation_time: "<<elapsed_time.diff<<" "<<elapsed_time.unit<<" ("<<elapsed_time.diff_fine<<" "<<elapsed_time.unit_fine<<")"<<'\n';
  log_file<<"mps: "<<moves_ps<<'\n';
  log_file<<"emps: "<<e_moves_ps<<'\n';
  log_file<<"fps: "<<flips_ps<<'\n';
  log_file<<"efps: "<<e_flips_ps<<'\n';
  log_file.close();
  return log_name;
}

};


#endif //NODES_HPP_FLIPPYBENCHMARKLOGGER_H
