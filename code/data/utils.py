import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
import re

import h5py
from tqdm import tqdm

def load_and_format_data_as_df(path, FPS=29.97):
    """load and format the hdf5 trajectory data as a single pandas dataframe with the following columns:
    time (float): time in seconds
    x (float): x position in meters
    y (float): y position in meters
    z (float): z position in meters
    vx (float): x velocity in meters per second
    vy (float): y velocity in meters per second
    vz (float): z velocity in meters per second
    V (float): velocity magnitude in meters per second
    tid (int): track identifier (unique for each particle in the entire dataset)

    Args:
        path (string): path to hdf5 data file
        FPS (float, optional): frame rate of the data recording equipment. Defaults to 29.97 Hz.

    Returns:
        _type_: _description_
    """
    with h5py.File(path, 'r') as file:
        data_all = {key: [] for key in file[list(file.keys())[0]]}
        data_all.update({'time': []})
        for frame_id in tqdm(file.keys()):
            time = np.ones(len(file[frame_id]['x'])) * float(frame_id) / FPS
            data_all['time'] = np.append(data_all['time'], time)
            for key in file[frame_id].keys():
                data_all[key] = np.append(data_all[key], file[frame_id][key])
        data_all = pd.DataFrame(data_all)
    return data_all

def get_birds_in_time_range(df, time_start, time_end):
    filtered_df = df[(df['time'] >= time_start) & (df['time'] <= time_end)].copy()
    times = filtered_df.time.unique()
    complete_tids = filtered_df.groupby('tid').filter(lambda x: x['time'].size == times.size)['tid'].unique()
    final_df = filtered_df[filtered_df['tid'].isin(complete_tids)]
    return final_df

def get_velocity_order_parameter(df):
    # calculates the the length of the average velocity vector for each time step
    num_birds = len(df.tid.unique())  # get number of birds for averaging
    df[['vx', 'vy', 'vz']] = df[['vx', 'vy', 'vz']].div(df['V'], axis=0)  # normalize the velocity vectors
    vel_order = df.groupby('time')[['vx', 'vy', 'vz']].sum().apply(np.linalg.norm, axis=1) / num_birds  # average the velocity vectors and get the resulting norm
    return vel_order