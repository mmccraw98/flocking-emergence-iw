%% start
clear; close all; clc;


L = 10;
num = 10;
Tmin = 50;
orders = zeros(num, 1);
etas = linspace(0.01, 3, num);
stds = zeros(num, 1);
means = zeros(num, 1);

for ii = 1:length(orders)
    [vel_order, time, x, y, vx, vy, std_v, mean_v] = vicsek_sim(L, etas(ii), 2000, 2, 1, 1, 100, false);
    orders(ii) = mean(vel_order(time > Tmin));
    stds(ii) = mean(std_v(time > Tmin));
    means(ii) = mean(mean_v(time > Tmin));
    etas(ii)
end


%%

% plot(etas, orders)
plot(etas, stds)


%%

[vel_order, time, x, y, vx, vy, std_v, mean_v] = vicsek_sim(L, etas(ii), 2000, 2, 1, 1, 100, false);



















%% functions

function [distance, dx, dy] = get_distances_pbc(x, y, L)
    dx = x - x';
    dy = y - y';
    
    % implement pbc
    dx = dx - L * round(dx / L);
    dy = dy - L * round(dy / L);
    
    distance = sqrt(dx .^ 2 + dy .^ 2);
end


function [vel_order, time, x, y, vx, vy, std_v, mean_v] = vicsek_sim(L, eta, N, v0, r_neigh, dt, TMAX, show_plot)

    neigh_avg = zeros(N, 1);
    
    x = rand(N, 1) * L;
    y = rand(N, 1) * L;
    theta = (2 * rand(N, 1) - 1) * pi;
    vel_order = zeros(TMAX, 1);
    time = (1:TMAX) * dt;
    std_v = zeros(TMAX, 1);
    mean_v = zeros(TMAX, 1);
    
    for t = 1:TMAX
        [distance, dx, dy] = get_distances_pbc(x, y, L);
        
        neigh_avg = neigh_avg * 0.0;
        for ii = 1:N
           neigh_dist = distance(ii, :);
           neigh_theta = theta(0 < neigh_dist & neigh_dist < r_neigh);
           if length(neigh_theta) > 1
               neigh_avg(ii) = atan2(mean(sin(neigh_theta)), mean(cos(neigh_theta)));
           else
               neigh_avg(ii) = theta(ii);
           end
        end
        
        vx = v0 * cos(theta);
        vy = v0 * sin(theta);
    
        x = x + vx * dt;
        y = y + vy * dt;
        % there is no need to update the positions as the periodic boundaries are
        % handled in the distance calculation
        
        theta = neigh_avg + (2 * rand(N, 1) - 1) * eta * pi / 2;
        % all the vicsek papers in literature that i have read use this
        % range of noise rather than the -1,1 which was given in the matlab
        % code linked in the homework - i will use the original values
    
        vel_order(t) = sqrt(mean(vx) ^ 2 + mean(vy) ^ 2) / mean(sqrt(vx .^ 2 + vy .^ 2));
        mean_v(t) = mean(sqrt((vx) .^ 2 + (vy) .^ 2));
        std_v(t) = std(sqrt((vx) .^ 2 + (vy) .^ 2));
        if show_plot
            quiver(mod(x, L), mod(y, L), vx, vy)
            xlim([0 L]);
            ylim([0 L]);
            axis square
            pause(0.1)
        end
    end
end
