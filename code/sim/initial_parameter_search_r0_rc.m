%% start
clear; close all; clc;
N = 100;  % num particles
dim_r0 = 2.0;  % ratio of line of sight to excluded volume (bird size r0 / rc)
dim_v0 = 1.0;  % dimensionless velocity magnitude (v0 dt / rc)
eta = 0.1;  % noise magnitude
beta = 5.0;  % excluded volume magnitude
theta = pi;  % line of sight arc
gamma = 0.1;  % strength of attraction to chimney
Nsteps = 1000;  % length of run
Nlog = 1000;  % number of steps to log
show_plot = true;  % show animation
save_last_N = 500;  % all the simulation data beyond this point will be saved

num_run = 10;
POL = zeros(num_run, 1);
MA = zeros(num_run, 1);
dim_r0s = linspace(1.01, 10.0, num_run);
dim_r0s(end + 1) = 20;
dim_r0s(end + 1) = 100;
dim_r0s(end + 1) = 500;
dim_r0s(end + 1) = 1000;
for jj = 1:length(dim_r0s)
    pol = zeros(Nsteps, 1);
    ma = zeros(Nsteps, 1);
    num = 10;
    dim_r0 = dim_r0s(jj);
    parfor ii = 1:num
        [xtotal, ytotal, vxtotal, vytotal, pos_cm, pols, mas] = vicsek_IW(N, dim_r0, dim_v0, eta, beta, gamma, theta, Nsteps, Nlog, show_plot, save_last_N);
        pol = pol + pols;
        ma = ma + mas;
    end
    pol = pol / num;
    ma = ma / num;
    POL(jj) = mean(pol(end - save_last_N: end));
    MA(jj) = mean(ma(end - save_last_N: end));
    jj, length(dim_r0s)
end
%%
clc; close;

figure (1); hold on;
plot(dim_r0s, POL, "Linewidth", 2);
plot(dim_r0s, MA, "Linewidth", 2);
legend("$\Phi$", "$\Psi$", 'Interpreter', 'Latex', 'Fontsize', 30);
xlabel("$r_0 / r_c$", 'Interpreter', 'Latex');
ylabel("Order Parameters", 'Interpreter', 'Latex');
set(gca, 'Xscale', 'Log', 'Fontsize', 20)
saveas(gcf, sprintf('ro_rc_analysis/eta%0.1f-beta%0.1f-gamma%0.1f-theta%0.1f.png', eta, beta, gamma, theta));

% under these settings, the ratio of r0/rc does not matter that much in the
% region where r0/rc > 1 - in the case of birds, it is unlikely that 
% r0/rc ~ 1, thus should be reasonable to say that r0/rc does not have an
% effect

name = sprintf('ro_rc_analysis/eta%0.1f-beta%0.1f-gamma%0.1f-theta%0.1f.csv', eta, beta, gamma, theta);
data = [POL MA dim_r0s'];
csvwrite(name, data);







