%% start
clear; close all; clc;

% parameters
dim_v0 = 1.03;  % dimensionless velocity magnitude
dim_eta = 1.0;  % dimensionless angular noise
phi = pi;  % rear field of view restriction
rho = 2.5;  % particle number density
% L = sqrt(N / rho);  % box length
% N = 2000;  % number of particles
L = 20;  % box length
N = L ^ 2 * rho;  % number of particles
N_steps = 2000;  % number of simulation steps
show_plot = true;  % whether to show movie while running

% fixed values
r = 1.0;  % field of view radius
dt = 1.0;  % default timestep
omega_max = 0.1745 / dt;  % maximum angular velocity (from main paper)
eta = dim_eta * omega_max * dt;  % noise strength
v0 = dim_v0 * omega_max * r;  % velocity magnitude



% positions
x = rand(1, N) * L;
y = rand(1, N) * L;

theta = rand(1, N) * 2 * pi;  % orientations
theta_neigh = zeros(1, N);  % neighbor orientations

ones_temp = ones(1, N);  % for vectorization

for tt = 1:N_steps
    % calculate distances using pbc convention
    [distance, dx, dy] = get_distances_pbc(x, y, L);
    
    % perform orientation averaging for all particles
    for ii = 1:N
        neigh_dist = distance(ii, :);
        
        dx_i = dx(ii, :);
        dy_i = dy(ii, :);
        
        % find if within field of view
        alpha = atan2(dy_i, dx_i);
        alpha = angdiff(ones_temp * theta(ii), alpha);
        index = (((neigh_dist <= r) & (abs(alpha) <= phi / 2)) | (neigh_dist == 0));
        
        % update neighbor velocities
        neigh_theta = theta(index);
        theta_neigh(ii) = atan2(mean(sin(neigh_theta)), mean(cos(neigh_theta)));
    end
    
    % calculate velocities
    vx = v0 * cos(theta);
    vy = v0 * sin(theta);
    
    % update positions
    x = x + vx * dt;
    y = y + vy * dt;
    
    % check if theta update (from neighbor average alone) exceeds speed limit
    dTheta = angdiff(theta, theta_neigh);
    theta_adj = (dTheta >= omega_max * dt) * omega_max * dt - (dTheta <= -omega_max * dt) * omega_max * dt;
    
    % update angles
    theta = theta_neigh + (2 * rand(1, N) - 1) * eta / 2;

    if anynan(theta)
        break;
    end
    if show_plot
        quiver(mod(x, L), mod(y, L), cos(theta), sin(theta));
        xlim([0 L]);
        ylim([0 L]);
        axis square;
        pause(0.01);
    end
end











