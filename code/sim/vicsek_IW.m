function [xtotal, ytotal, vxtotal, vytotal, pols, mas] = vicsek_IW(N, pf, r0, v0, dt, eta, beta, phi_los, gamma, Nsteps, Nplot, show_plot, save_last_N)

% Box Size
L = 20.0;

% Definte particle sizes (rc) based on packing fraction
rc = 2.0 * L * sqrt(pf / (pi * N));

% scale r0 by rc
r0 = r0 * rc;

% initialize trajectory values
xtotal = zeros(N, save_last_N);
ytotal = zeros(N, save_last_N);
vxtotal = zeros(N, save_last_N);
vytotal = zeros(N, save_last_N);

% initialize order parameter values
pols = zeros(Nsteps, 1);
mas = zeros(Nsteps, 1);

% Initialize Positions and Velocities randomly
if gamma == 0
    rs = rand(N, 2) * L;
else
    rs = (rand(N, 2) - 0.5) * 2 * L;
end
vs = randn(N, 2);
vnorm = sqrt(sum(vs'.^2))';
vs = vs .* v0 ./ [vnorm, vnorm];

% Loop over simulation
for step = 1:Nsteps
    % Perform Integration
    vs = vicsekvelocities_IW(N, v0, r0, rc, eta, beta, L, rs, vs, phi_los, gamma) ;
    rs = rs + vs * dt;
    
    % store order parameters
    meanVs = mean(vs);
    meanVsNrm = sqrt(meanVs(1) ^ 2 + meanVs(2) ^ 2);
    pols(step) = meanVsNrm / v0;
    rcm = rs - mean(rs, 1);
    mas(step) = mean(abs(rcm(:, 1) .* vs(:, 2) / v0 - rcm(:, 2) .* vs(:, 1) / v0) ./ sqrt(sum(rcm' .^ 2))');
    
    if step > Nsteps - save_last_N
        % store trajectories and velocities
        xtotal(:,step) = rs(:,1);
        ytotal(:,step) = rs(:,2);
        vxtotal(:,step) = vs(:,1);
        vytotal(:,step) = vs(:,2);
    end


    
    % Plot cells and velocities
    if mod(step, Nplot) == 0
        fprintf('On step %d for beta = %f, eta = %f\n',step,beta,eta);
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