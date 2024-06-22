import os
import sys
from typing import List, Dict, NamedTuple
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9
import yaml

LogPathData = namedtuple('LogData', 'path timestamp')

ParsedData = Dict[str, str]


class LogData(NamedTuple):
    timestamp: str
    parsed_data: ParsedData


def get_subdirectories(directory: str) -> ([str], [str]):
    sub_dir_paths = []
    sub_dir_names = []
    for d in os.listdir(directory):
        full_path = os.path.join(directory, d)
        if os.path.isdir(full_path):
            sub_dir_paths.append(full_path)
            sub_dir_names.append(d)
    return sub_dir_paths, sub_dir_names


def get_files(directory: str, day: str) -> List[LogPathData]:
    path_data: List[LogPathData] = []
    for f in os.listdir(directory):
        full_path = os.path.join(directory, f)
        if os.path.isfile(full_path) and f.endswith('.yml'):
            hour = f.split('.')[0]
            data = LogPathData(path=full_path, timestamp=f'{day}-{hour}')
            path_data.append(data)
    return path_data


def _parse_log(log: str)-> Dict[str, str]:
    return yaml.safe_load(log)


def parse_log(log_path_data: LogPathData):
    log_data: LogData
    with open(file=log_path_data.path, mode='r') as log_file:
        raw_data = log_file.read()
        log_data = LogData(parsed_data=_parse_log(raw_data),
                           timestamp=log_path_data.timestamp)
    return log_data


def get_all_logs(logs_dir: str,
                 benchmark_version: int | None = None) -> List[ParsedData]:

    assert os.path.isdir(logs_dir), f'{logs_dir=} is not a valid path'

    log_sub_dirs, log_dates = get_subdirectories(logs_dir)

    all_log_paths: List[LogPathData] = []
    for lsdp, day in zip(log_sub_dirs, log_dates):
        files = get_files(lsdp, day)
        all_log_paths += files

    all_logs: List[ParsedData] = []
    for log_path_data in all_log_paths:
        log = parse_log(log_path_data=log_path_data)
        if benchmark_version is not None:
            if ("BENCHMARK VERSION" not in log.parsed_data or
                    log.parsed_data['BENCHMARK VERSION'] != benchmark_version):
                continue
        log.parsed_data['timestamp'] = log.timestamp
        all_logs.append(log.parsed_data)
    return all_logs


def main():
    logs_dir = sys.argv[1]
    all_logs: List[ParsedData] = get_all_logs(logs_dir)

    df = pd.DataFrame(all_logs)
    def version_num_to_category(num: float):
        if not pd.isna(num):
            return str(int(num))
        else:
            return '-'

    df['BENCHMARK VERSION'] = df['BENCHMARK VERSION'].apply(version_num_to_category)
    df['timestamp'] = pd.to_datetime(df['timestamp'], format='%Y-%m-%d-%H-%M-%S')
    df.sort_values(by='timestamp')


    df_long = pd.melt(df, id_vars=['timestamp', 'BENCHMARK VERSION'],
                      value_vars=['mps', 'emps', 'fps', 'efps'],
                      var_name='perf_type', value_name='perf_value')
    df_long.sort_values(by='timestamp')
    df_long['timestamp'] = df_long['timestamp'].dt.strftime('%Y-%m-%d %H:%M')
    print(df_long.tail())

    # x_breaks = df_long['timestamp'][::2].to_list()
    # print(f'{x_breaks=}')
    plot = (p9.ggplot(df_long)
            + p9.aes(x='timestamp', y='perf_value', color='BENCHMARK VERSION')
            # + p9.geom_text(p9.aes(label='BENCHMARK VERSION'), size=9)
            + p9.geom_point()
            + p9.facet_wrap('~perf_type', scales='free_y')
            + p9.theme_minimal()
            + p9.labs(x='time', y='perf', color='Bench.\nVersion',
                      title='performance over time')
            # + p9.scale_color_discrete(name='Bench.\nV')
            # + p9.scale_x_discrete(breaks=x_breaks)
            # + p9.scale_x_datetime(breaks=x_breaks)
            + p9.theme(
                figure_size=(15, 9),
                axis_text_x=p9.element_text(angle=70, hjust=1),
            )
            )

    print(plot)


if __name__ == '__main__':
    main()