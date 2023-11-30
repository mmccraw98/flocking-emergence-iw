function [distance, dx, dy] = get_distances_pbc(x, y, L)
    dx = x - x';
    dy = y - y';
    
    % implement pbc
    dx = dx - L * round(dx / L);
    dy = dy - L * round(dy / L);
    
    distance = sqrt(dx .^ 2 + dy .^ 2);
end