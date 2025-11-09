function [area_layer, spad_layer, spad_mesh] = distribute_spad_mesh_by_layer_decay( ...
        group_sp_t, group_se, scale, layer_thickness, total_spad, Kn)
%% DISTRIBUTE_SPAD_MESH_BY_LAYER_DECAY
% Distribute total leaf SPAD to each mesh (triangular facet) in a canopy based on vertical layers
% and an exponential decay weighting.
%
% Inputs:
%   group_sp_t      V×3 vertex coordinates of the canopy
%   group_se        F×3 face indices (triangles)
%   scale           scaling factor (consistent with leaf area calculation)
%   layer_thickness thickness of each vertical layer
%   total_spad      total SPAD value to distribute
%   Kn              decay coefficient (for vertical distribution)
%
% Outputs:
%   area_layer      1×n_layers, total leaf area of each layer
%   spad_layer      1×n_layers, SPAD value assigned to each layer (sum = total_spad)
%   spad_mesh       F×1, SPAD value assigned to each triangular mesh (sum = total_spad)

    %% 1. Scale vertex coordinates and compute area of each triangle
    p_scaled   = group_sp_t .* scale ./ 100;           % scale coordinates
    leaf_area  = LA(group_sp_t, group_se, scale);     % compute leaf area for each triangle

    %% 2. Determine number of vertical layers & actual layer thickness
    z_all      = p_scaled(:,3);                        % all z-coordinates
    z_min      = min(z_all);
    z_max      = max(z_all);
    n_layers   = ceil((z_max - z_min) / layer_thickness); % number of vertical layers
    dz         = (z_max - z_min) / n_layers;          % adjusted layer thickness

    %% 3. Compute centroid z-coordinate of each triangle & assign to layer
    z1         = p_scaled(group_se(:,1),3);
    z2         = p_scaled(group_se(:,2),3);
    z3         = p_scaled(group_se(:,3),3);
    z_centroid = (z1 + z2 + z3) / 3;                 % triangle centroid
    layer_raw  = floor((z_centroid - z_min) ./ dz) + 1; % raw layer index
    layer_raw  = max(1, min(n_layers, layer_raw));     % ensure within bounds

    %% 4. Flip layer index so top layer = 1, bottom layer = n_layers
    layer_idx  = n_layers - layer_raw + 1; 

    %% 5. Sum leaf area per layer
    area_layer = accumarray(layer_idx, leaf_area, [n_layers, 1])';  % 1×n_layers

    %% 6. Compute normalized layer depth x ∈ (0,1) and exponential weight N_pct
    x_layer    = ((1:n_layers) - 0.5) ./ n_layers;       % normalized depth of layer center
    N_pct_layer = Kn .* exp(-Kn .* x_layer);             % weight for each layer

    %% 7. Compute weighted SPAD for each layer
    weighted_area = area_layer .* N_pct_layer;           
    spad_layer    = total_spad .* (weighted_area / sum(weighted_area)); % SPAD per layer

    %% 8. Distribute layer SPAD to individual triangles proportionally to area
    spad_density  = spad_layer ./ area_layer;           % SPAD per unit area
    spad_mesh     = spad_density(layer_idx)' .* leaf_area;  % SPAD per triangle

    % Optional validation
    % assert( abs(sum(spad_layer)-total_spad) < 1e-6 );
    % assert( abs(sum(spad_mesh)-total_spad)  < 1e-6 );
    abs(sum(spad_mesh)-total_spad)
end


%% Auxiliary function: compute triangle areas
function leaf_area = LA(vertices, faces, BL)
    p1 = vertices(faces(:,1), :);
    p2 = vertices(faces(:,2), :);
    p3 = vertices(faces(:,3), :);

    d1 = sqrt(sum((p1 - p2).^2, 2));
    d2 = sqrt(sum((p1 - p3).^2, 2));
    d3 = sqrt(sum((p3 - p2).^2, 2));

    % Heron's formula
    pp = (d1 + d2 + d3) / 2;
    leaf_area = sqrt(abs(pp .* (pp - d1) .* (pp - d2) .* (pp - d3))) * BL^2;
end
