mkdir macro_bench_build
cd macro_bench_build || exit

readonly bench_dir="../benchmarks/rbc_bench"
cmake ${bench_dir}

make

jq '.save_dir="../demo_out/rbc_bench" | .log_dir="../logs/rbc_bench"' "${bench_dir}/rbc_config.json" > rbc_config.json

./flippy_rbc_benchmark rbc_config.json

rm macro_bench_build
