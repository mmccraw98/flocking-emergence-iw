%% start
clear; close all; clc;
N = 1800;  % num particles
pf = 0.2;  % packing fraction
r0 = 2.0;  % ratio of line of sight to excluded volume (bird size)
v0 = 1.0;  % velocity magnitude
dt = 1.0;  % timestep
eta = 0.5;  % noise magnitude
beta = 0.2;  % excluded volume magnitude
phi_los = pi;  % line of sight cone
gamma = 0.2;  % strength of attraction to chimney
Nsteps = 1000;  % length of run
Nplot = 100;  % number of steps to log
show_plot = true;  % show animation
save_last_N = 500;  % all the simulation data beyond this point will be saved
[xtotal, ytotal, vxtotal, vytotal, pols, mas] = vicsek_IW(N, pf, r0, v0, dt, eta, beta, phi_los, gamma, Nsteps, Nplot, show_plot, save_last_N);

%%
clc;

figure (1); hold on;
plot(pols)
plot(mas)
legend('Velocity Order', 'Angular Momentum Order')
saveas(gcf, 'sim_order_sample.png')













