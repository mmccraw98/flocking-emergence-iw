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

def get_principal_axes(df_start, skip_every_n_steps=3):
    pca = PCA(n_components=3)

    principal_axes = []

    df = df_start.copy()
    times = np.sort(df.time.unique())
    for i in tqdm(range(1, times.size // skip_every_n_steps)):
        time = times[i * skip_every_n_steps]
        sub_df = df[(df.time == time)].copy()
        pos = sub_df[['x', 'y', 'z']].values
        try:
            pca.fit(pos)
            principal_axes.append(pca.components_)
        except:
            principal_axes.append(np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]))
    
    return np.array(principal_axes)

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
        if sdf.shape[0] < 10:
            continue
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
        points = df[df.time == df.time.values[i]][['x', 'y', 'z']].values
        dists = get_dists(points)
        bin_values += np.histogram(dists, bins=bins)[0]

    pdf = bin_values / (np.sum(bin_values) * (bin_centers[1] - bin_centers[0]))
    return bin_centers, pdf

def np_convolve_moving_average(arr, window_size):
    return np.convolve(arr, np.ones(window_size) / window_size, mode='valid')

def section_df_into_disordered_ordered_regimes(df_start, phi_low=0.2, phi_high=0.6, sm_perc=0.01, min_length=100):
    df = df_start.copy()
    vel_order = get_velocity_order_parameter(df)
    time = np.sort(df.time.unique())

    win_size = int(vel_order.size * sm_perc)
    vel_order_sm = np_convolve_moving_average(vel_order, win_size)

    results = {'low': {'regimes': [], 'condition': vel_order_sm < phi_low}, 'high': {'regimes': [], 'condition': vel_order_sm > phi_high}}

    for key in results.keys():
        condition = results[key]['condition']
        starts = np.where(np.diff(condition.astype(int)) == 1)[0] + 1
        stops = np.where(np.diff(condition.astype(int)) == -1)[0] + 1

        if condition[0]:
            starts = np.insert(starts, 0, 0)
        if condition[-1]:
            stops = np.append(stops, len(condition))

        offset = int(win_size / 2)
        adjusted_starts = starts + offset
        adjusted_stops = stops + offset

        adjusted_starts = np.clip(adjusted_starts, 0, len(vel_order) - 1)
        adjusted_stops = np.clip(adjusted_stops, 0, len(vel_order))

        for start, stop in zip(adjusted_starts, adjusted_stops):
            if stop - start > min_length:
                results[key]['regimes'].append(time[start:stop])
    return results['low']['regimes'], results['high']['regimes']

def time_correlation_function(vectors, norm=True):
    N = len(vectors)
    max_lag = N - 1
    correlation = np.zeros(max_lag)

    # Compute the autocorrelation for delta t = 0
    if norm:
        autocorr_at_zero = np.mean([np.dot(vectors[t], vectors[t]) for t in range(N)])
    else:
        autocorr_at_zero = 1

    for dt in tqdm(range(max_lag)):
        sum_corr = 0
        count = 0
        for t in range(N - dt):
            sum_corr += np.dot(vectors[t], vectors[t + dt])
            count += 1
        correlation[dt] = sum_corr / (count * autocorr_at_zero)

    return correlation
