#! python3
import sys
import json
import re
import numpy as np
from termcolor import colored, cprint

if __name__ == '__main__':
    bench_types = dict()
    path = sys.argv[1]
    bench_file = sys.argv[2]
    with open(f"{path}/{bench_file}", 'r') as f:
        j = json.load(f)
    reg = re.compile(r'(^[a-zA-Z][\w\s]*)(base|experimental)', re.IGNORECASE)

    for benchmark in j['benchmarks']:
        full_name = benchmark["name"]
        name = reg.match(full_name).group(1)
        suffix = reg.match(full_name).group(2)
        if not bench_types.get(name):
            bench_types[name] = dict()
        if benchmark['aggregate_name'] == 'mean':
            bench_types[name][suffix] = {"mean": benchmark['real_time']}
        if benchmark['aggregate_name'] == 'stddev':
            bench_types[name][suffix]["std_dev"] = benchmark['real_time']

    diff_log_name = "bench_diff.log"
    diff_log_file = open(f"{path}/bench_diff.log", "w")
    diff_json_file = open(f"{path}/bench_diff.json", "w")

    diff_dict = dict()

    string_stencil = f"{'test name':>30}{' ':<3}{'base time':<9}  |  {'expr time':<9}  |  {'rel impr':<9}\n"
    line_width=len(string_stencil)
    string_stencil += f"{' ':>30}{' ':<3}{'base std':<9}  |  {'expr std':<9}  |  {'av_error/tot_diff':<9}\n"
    string_stencil += line_width*'_'+'\n'
    for benchmark_name, benchmark in bench_types.items():
        base_time = benchmark['base']['mean']
        base_std = benchmark['base']['std_dev']
        experimental_time = benchmark['experimental']['mean']
        experimental_std = benchmark['base']['std_dev']
        diff = base_time - experimental_time
        av_noise = (base_std + experimental_std)/2
        rel_change = diff/base_time
        rel_noise = av_noise / abs(diff)

        if not diff_dict.get(benchmark_name):
            diff_dict[benchmark_name] = dict()

        diff_dict[benchmark_name]["base time"] = base_time
        diff_dict[benchmark_name]["base std"] = base_std

        diff_dict[benchmark_name]["experimental time"] = experimental_time
        diff_dict[benchmark_name]["experimental std"] = experimental_std

        diff_dict[benchmark_name]["diff"] = diff
        diff_dict[benchmark_name]["av_noise"] = av_noise

        diff_dict[benchmark_name]["rel change"] = rel_change
        diff_dict[benchmark_name]["rel error"] = rel_noise

        attr = ['bold', 'dark']
        if rel_noise < 1:
            if diff > 0:
                ticker = '∆∆∆'
                color = 'cyan'
                if rel_noise<0.1:
                    color = 'green'
                    attr.append('blink')
            else:
                ticker = 'VVV'
                color = 'red'
                attr.append('blink')
        else:
            ticker = '---'
            color = 'yellow'

        res = f'{np.round(rel_change,3)} {colored(ticker, color, attrs=attr)}'
        string_stencil += f"{benchmark_name:>30}{' ':<3}{np.round(base_time,3):<9}  |  {np.round(experimental_time,3):<9}  |  {res:<9}\n"
        string_stencil += f"{' ':>30}{' ':<3}{np.round(base_std,3):<9}  |  {np.round(experimental_std,3):<9}  |  {np.round(rel_noise,3):<9}\n"
        string_stencil += int(line_width/4)*' -- '+'\n'

    print(string_stencil, file=diff_log_file)
    print(string_stencil)

    out_bench = {'context': j['context'], 'diffs': diff_dict }
    json.dump(out_bench, diff_json_file, indent = 6)
