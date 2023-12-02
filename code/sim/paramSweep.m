%% start
clear; close all; clc;
N = 1000;  % num particles
Nsteps = 1000;  % length of run
Nlog = 1001;  % number of steps to log
show_plot = false;  % show animation
save_last_N = 500;  % all the simulation data beyond this point will be saved

theta = 2 * pi;  % line of sight arc - ADJUST BETWEEN PI AND 2 PI (for vicsek)

dim_r0s = [2 5 10];
dim_v0s = [1e-1 1e0 1e1];
log_axis = [1e-3 1e-2 1e-1 1e0];  % for eta, beta, gamma % MAYBE ADD ANOTHER AT 1e0!

num_repeats = 10;  % each simulation will be repeated this many times and the results will be averaged

count = 1;
total = length(dim_r0s) * length(dim_v0s) * length(log_axis) ^ 3;
for i1 = 1:length(dim_r0s)
    for i2 = 1:length(dim_v0s)
        for i3 = 1:length(log_axis)  % eta
            for i4 = 1:length(log_axis)  % beta
                for i5 = 1:length(log_axis)  % gamma
                    % unpack the parameters now to avoid parfor overhead
                    dim_r0 = dim_r0s(i1);  % ratio of line of sight to excluded volume (bird size r0 / rc)
                    dim_v0 = dim_v0s(i2);  % dimensionless velocity magnitude (v0 dt / rc)
                    eta = log_axis(i3);  % noise magnitude
                    beta = log_axis(i4);  % excluded volume magnitude
                    gamma = log_axis(i5);  % strength of attraction to chimney

                    % check if the path exists to avoid overwriting
                    ordername = sprintf('param_sweep/eta1e%d_beta1e%d_gamma1e%d_theta%dpi_ro%d_vo1e%d_.csv', log10(eta), log10(beta), log10(gamma), theta / pi, dim_r0, log10(dim_v0));
                    if isfile(ordername)
                        fprintf(['Skipping ' ordername '\n']);
                        count = count + 1;
                        continue
                    end
                    tic;

                    
                    % do the ensemble averaging
                    pol = zeros(Nsteps, 1);
                    ma = zeros(Nsteps, 1);
                    parfor ii = 1:num_repeats
                        [xtotal, ytotal, vxtotal, vytotal, pos_cm, pols, mas] = vicsek_IW(N, dim_r0, dim_v0, eta, beta, gamma, theta, Nsteps, Nlog, show_plot, save_last_N);
                        pause(0.01 * (ii - 1));  % to avoid parallel overwrite - maybe...
                        pol = pol + pols;
                        ma = ma + mas;
                        
                        % save the center of mass data
                        name = sprintf('param_sweep_com_pos/eta1e%d_beta1e%d_gamma1e%d_theta%dpi_ro%d_vo1e%d_run%d_.csv', log10(eta), log10(beta), log10(gamma), theta / pi, dim_r0, log10(dim_v0), ii);
                        data_as_table = array2table(pos_cm');
                        data_as_table.Properties.VariableNames(1:2) = {'xcm', 'ycm'};
                        writetable(data_as_table, name);
                    end
                    pol = pol / num_repeats;
                    ma = ma / num_repeats;

                    % save the results
                    data = [pol ma];
                    data_as_table = array2table(data);
                    data_as_table.Properties.VariableNames(1:2) = {'va', 'ma'};
                    writetable(data_as_table, ordername);

                    fprintf('__________________________________\nDone with %d of %d\nIteration time: %0.1f s\nRemaining time: %f s\n', count, total, toc, toc * (total - count));
                    count = count + 1;
                end
            end
        end
    end
end

