function vs = vicsekvelocities_IW(N, v0, r0, rc, eta, beta, gamma, theta, rs, vs)
    %{
    Function that finds the updated velocity for a Vicsek model with
    repulsion using the scheme from PhysRevE.77.046113 with a modification
    to allow line of sight alignment from 10.1088/1361-6463/aab0d4 and a
    global attraction towards the chimney from PhysRevE.107.014209
    
    Inputs:
        N          The number of cells
        v0         The constant speed of the cells
        r0         The distance cutoff within which cells are considered neighbors 
        rc         The diameter of the cell
        eta        The relative strength of the noise
        beta       The relative strength of the repulsion
        gamma      The strength of the attraction to the origin
        theta      The angle of the particle line of sight for neighbor alignment
        rs         Nx2 matrix of cell positions
        vs         Nx2 matrix of cell velocities
    
    Outputs:
        vs     Nx2 matrix of updated velocities due to the repulsive Vicsek algorithm
    
    %}
    
    % Initialize Variables
    sum_vs = zeros(N, 2);   % Sum of neighboring velocities
    Si_norm = zeros(N, 1);  % Number of neighbors for each cell i
    Fi = zeros(N, 2);       % Total repulsive force for each cell i
    ones_temp = ones(N, 1);
    orientation = atan2(vs(:, 2), vs(:, 1));

    % Calculate Fi
    for i=1:N
       % Calculate distance between cell i and the rest within a periodic box
       rijs = -[rs(:, 1) - rs(i, 1), rs(:, 2) - rs(i, 2)];

       dists = sqrt(sum(rijs' .^ 2))';

       % calculate the norm for direction
       rijs_norm = normer(rijs);
       
       % Calculate the set Si using the distance and angle constraint
       alpha = atan2(rijs_norm(:, 2), rijs_norm(:, 1));
       alpha = angdiff(ones_temp * orientation(i), alpha);
       Si = (((dists <= r0) & (abs(alpha) <= theta / 2)) | (dists == 0));
       
       Si_norm(i) = sum(Si);
       
       % Calculate sum_vs(i,:)
       sum_vs(i, :) = sum(vs(Si, :), 1);
       
       % Calculate the repulsive force due to each cell on i
       index = (dists <= rc) & (dists > 0);
       force_mag = (1 - dists(index) / rc);
       
       Fi(i,:) = sum(rijs_norm(index, :) .* force_mag, 1);
    end
    
    % Calculate a matrix of random unit vectors
    noise = normer(randn(N, 2));
    
    vs = normer(sum_vs ./ v0 + beta .* Fi + eta .* [Si_norm, Si_norm] .* noise - gamma * rs) .* v0;
end


function v = normer(input)
    %{
    A function that divides each row in a matrix by its norm.
    
    Inputs:
        input     A matrix for which to normalize along the second axis
    Outputs:
        v         The normalized matrix
    %}
    
    normval = sqrt(sum(input .^ 2, 2)) + 10 ^ - 16;
    v = input ./ [normval, normval];
end