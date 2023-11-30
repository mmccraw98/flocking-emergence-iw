%% start
clear; close all; clc;
N = 1800;  % num particles
pf = 0.5;  % packing fraction
r0 = 2.0;  % ratio of line of sight to excluded volume (bird size)
v0 = 1.0;  % velocity magnitude
dt = 1.0;  % timestep
eta = 0.2;  % noise magnitude
beta = 0.5;  % excluded volume magnitude
phi_los = pi;  % line of sight cone
gamma = 1.0;  % strength of attraction to chimney
Nsteps = 2000;  % length of run
Nplot = 100;  % number of steps to log
show_plot = true;  % show animation
[xtotal, ytotal, pols] = vicsek_IW(1800,0.7,2.0,0.2,1.0,0.2,0.5, pi, 1.0, 2000,100,true)