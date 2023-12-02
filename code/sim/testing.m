%% start
clear; close all; clc;
N = 1000;  % num particles
Nsteps = 1000;  % length of run
Nlog = 1001;  % number of steps to log
show_plot = true;  % show animation
save_last_N = 500;  % all the simulation data beyond this point will be saved

theta = pi;  % line of sight arc - ADJUST BETWEEN PI AND 2 PI (for vicsek)

% mill high, vel low
gamma = 1e-2;
beta = 1e-3;
dim_v0 = 1e1;
dim_r0 = 2.0;
eta = 1e-3;

% mill med, vel high
gamma = 1e0;
beta = 1e0;
dim_v0 = 1e0;
dim_r0 = 2.0;
eta = 1e-1;


% mill low, vel high
gamma = 1e0;
beta = 1e-3;
dim_v0 = 1e-1;
dim_r0 = 2.0;
eta = 1e-3;

[xtotal, ytotal, vxtotal, vytotal, pos_cm, pols, mas] = vicsek_IW(N, dim_r0, dim_v0, eta, beta, gamma, theta, Nsteps, Nlog, show_plot, save_last_N);


                        

%%

plot(mas)







