function [xtotal, ytotal, vxtotal, vytotal, pos_cm, pols, mas] = vicsek_IW(N, dim_r0, dim_v0, eta, beta, gamma, theta, Nsteps, Nlog, show_plot, save_last_N)
%{
Function that simulates a modified Vicsek model with
repulsion using the scheme from PhysRevE.77.046113 with a modification
to allow line of sight alignment from 10.1088/1361-6463/aab0d4 and a
global attraction towards the chimney from PhysRevE.107.014209

Inputs:
    N          The number of cells
    dim_r0     The dimensionless line of sight radius (ratio of r0 / rc)
    dim_v0     The dimensionless velocity magnitude (ratio of v0 dt / rc)
    eta        The relative strength of the noise
    beta       The relative strength of the repulsion
    gamma      The strength of the attraction to the origin
    theta      The angle of the particle line of sight for neighbor alignment
    rs         Nx2 matrix of cell positions
    vs         Nx2 matrix of cell velocities
    Nsteps     Number of steps to simulate
    Nlog       Output a console message every Nlog steps
    save_N     Save the configuration for only the last N steps

Outputs:
    vs     Nx2 matrix of updated velocities due to the repulsive Vicsek algorithm

%}

% fixed values gives dimension
rc = 1.0;  % repulsion radius
dt = 1.0;  % time step
packing_fraction_initial = 0.5;  % this might be a pointless initial condition

% give scale to dimensionless variables
r0 = dim_r0 * rc;
v0 = dim_v0 * rc / dt;
L = sqrt(N / packing_fraction_initial);


% initialize trajectory values
xtotal = zeros(N, save_last_N);
ytotal = zeros(N, save_last_N);
vxtotal = zeros(N, save_last_N);
vytotal = zeros(N, save_last_N);
pos_cm = zeros(2, save_last_N);

% initialize order parameter values
pols = zeros(Nsteps, 1);
mas = zeros(Nsteps, 1);

% Initialize Positions and Velocities randomly
rs = (rand(N, 2) - 0.5) * 2 * L;

vs = randn(N, 2);
vnorm = sqrt(sum(vs'.^2))';
vs = vs .* v0 ./ [vnorm, vnorm];

% Loop over simulation
for step = 1:Nsteps
    % Perform Integration
    vs = vicsekvelocities_IW(N, v0, r0, rc, eta, beta, gamma, theta, rs, vs);
    rs = rs + vs * dt;
    
    % store order parameters
    meanVs = mean(vs);
    meanVsNrm = sqrt(meanVs(1) ^ 2 + meanVs(2) ^ 2);
    pols(step) = meanVsNrm / v0;
    rcm = rs - mean(rs, 1);
    mas(step) = mean(abs(rcm(:, 1) .* vs(:, 2) / v0 - rcm(:, 2) .* vs(:, 1) / v0) ./ sqrt(sum(rcm' .^ 2))');
    
    if step > Nsteps - save_last_N
        % store trajectories and velocities
        save_step = abs(Nsteps - step - save_last_N);
        xtotal(:, save_step) = rs(:, 1);
        ytotal(:, save_step) = rs(:, 2);
        vxtotal(:, save_step) = vs(:, 1);
        vytotal(:, save_step) = vs(:, 2);
        pos_cm(:, save_step) = mean(rs, 1);
    end


    
    % Plot cells and velocities
    if mod(step, Nlog) == 0
        fprintf('On step %d of %d', step, Nsteps);
    end

    if show_plot
        if gamma == 0        
            quiver(mod(rs(:, 1), L), mod(rs(:, 2), L), vs(:, 1), vs(:, 2));
            xlim([0 L]);
            ylim([0 L]);
        else
            quiver(rs(:, 1), rs(:, 2), vs(:, 1), vs(:, 2));
            xlim([min(-L, min(rs(:, 1))) max(L, max(rs(:, 1)))])
            ylim([min(-L, min(rs(:, 2))) max(L, max(rs(:, 2)))])
        end
        axis square;
        pause(0.01);
    end
end