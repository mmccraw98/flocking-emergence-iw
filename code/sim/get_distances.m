function [distance, dx, dy] = get_distances(x, y)
    dx = x - x';
    dy = y - y';
    
    distance = sqrt(dx .^ 2 + dy .^ 2);
end