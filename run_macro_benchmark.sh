mkdir macro_bench_build
cd macro_bench_build || exit

readonly bench_dir="../benchmarks/rbc_bench"
cmake ${bench_dir}

make

jq '.save_dir="../demo_out/rbc_bench/" | .log_dir="../logs/rbc_bench/"' "${bench_dir}/rbc_config.json" > rbc_config.json

./flippy_rbc_benchmark rbc_config.json

cd .. || exit
rm -r macro_bench_build

cd "demo_out/rbc_bench" || exit

readonly viz_py_folder="../../demo/biconcave_shapes_MC/"
readonly viz_py="data_visualisation.py"


cp "${viz_py_folder}/${viz_py}" ${viz_py}

python3 ${viz_py}

cd ../../benchmarks/ || exit 1

python3 make_benchmark_timeline.py ../logs/rbc_bench