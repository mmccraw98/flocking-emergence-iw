import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d

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

def get_velocity_order_parameter(df_start):
    df = df_start.copy()
    # calculates the the length of the average velocity vector for each time step
    num_birds = df.groupby('time').apply(lambda x: x.tid.nunique())  # get number of birds for averaging
    df[['vx', 'vy', 'vz']] = df[['vx', 'vy', 'vz']].div(df['V'], axis=0)  # normalize the velocity vectors
    vel_order = df.groupby('time')[['vx', 'vy', 'vz']].sum().apply(np.linalg.norm, axis=1) / num_birds  # average the velocity vectors and get the resulting norm
    return vel_order

def get_variance_ratios(df_start, skip_every_n_steps=3):
    pca = PCA(n_components=3)

    variance_ratio_components = []

    df = df_start.copy()
    times = np.sort(df.time.unique())
    for i in tqdm(range(1, times.size // skip_every_n_steps)):
        time = times[i * skip_every_n_steps]
        sub_df = df[(df.time == time)].copy()
        pos = sub_df[['x', 'y', 'z']].values
        try:
            pca.fit(pos)
            variance_ratios = pca.explained_variance_ratio_  # if the third component is large, then the birds are flying in a volume, if it is small, they are in a plane
            variance_ratio_components.append(variance_ratios)
        except:
            variance_ratio_components.append([0, 0, 0])  # sometimes there is no data, probably need to handle this better
    
    return np.array(variance_ratio_components)

def get_number_densities(df_start, skip_every_n_steps=3):
    pca = PCA(n_components=3)

    densities = []  # transform data to PCA space, get volume of the bounding ellipsoid, and divide by number of birds

    df = df_start.copy()
    times = np.sort(df.time.unique())
    for i in tqdm(range(1, times.size // skip_every_n_steps)):
        time = times[i * skip_every_n_steps]
        sub_df = df[(df.time == time)].copy()
        pos = sub_df[['x', 'y', 'z']].values
        try:
            pca_pos = pca.fit_transform(pos)
            ranges = np.ptp(pca_pos, axis=0)
            a, b, c = ranges / 2  # Half the lengths of principal axes
            volume = 4/3 * np.pi * a * b * c  # treat the flock as an ellipsoid
            densities.append(sub_df.tid.nunique() / volume)
        except:
            densities.append(0)  # sometimes there is no data, probably need to handle this better
    
    return np.array(densities)

def get_number_densities_alt(df_start, skip_every_n_steps=3):
    df = df_start.copy()
    densities = []
    times = np.sort(df.time.unique())
    for i in tqdm(range(1, times.size // skip_every_n_steps)):
        time = times[i * skip_every_n_steps]
        sub_df = df[(df.time == time)].copy()
        num_birds = sub_df.tid.nunique()
        try:
            cuh = ConvexHull(sub_df[['x', 'y', 'z']].values)
            densities.append(num_birds / cuh.volume)
        except:
            densities.append(0)
    return np.array(densities)

def downsample_A2B(A, B):
    x = np.linspace(0, 1, len(A))
    x_new = np.linspace(0, 1, len(B))
    f = interp1d(x, A)
    return f(x_new)

def get_dists(points):
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=-1))
    return distances[np.triu_indices_from(distances, k=1)]

def random_sample_upper_bounds(df, n_random=10):
    max_dist = -np.inf
    for i in range(n_random):
        sdf = df[df.time == df.time.unique()[np.random.randint(0, df.time.unique().size)]].copy()
        points = sdf[['x', 'y', 'z']].values
        dists = get_dists(points)
        max_dist = max(max_dist, dists.max())
    return max_dist

def get_dist_pdf(df, N_bins=100):
    sdf = df.copy()
    max_dist = random_sample_upper_bounds(df, n_random=10)
    bins = np.linspace(0, max_dist, N_bins)
    bin_centers = (bins[1:] + bins[:-1]) / 2
    bin_values = np.zeros(N_bins - 1)

    for i in tqdm(range(sdf.time.unique().size)):
        points = sdf[['x', 'y', 'z']].values
        dists = get_dists(points)
        bin_values += np.histogram(dists, bins=bins)[0]

    pdf = bin_values / (np.sum(bin_values) * (bin_centers[1] - bin_centers[0]))
    return bin_centers, pdf
