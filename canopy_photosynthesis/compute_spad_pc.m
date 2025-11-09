function [leaf_spad_mesh, plant_face] = compute_spad_pc(pc, tri, a, b)
%% ========================================================================
% Function: compute_spad_pc
% Purpose : Compute SPAD values for each mesh facet based on DGCI (Dark Green Color Index)
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
% Inputs  :
%    pc  - pointCloud object with color information
%    tri - Nx3 array of mesh faces (indices into pc.Location)
%    a,b - linear coefficients for SPAD = a * DGCI - b
% Outputs :
%    leaf_spad_mesh - median SPAD value for each mesh face
%    plant_face     - filtered faces after removing noisy points
%% ========================================================================

% Extract color data
color = double(pc.Color);  
num_points = length(color);

% Pre-allocate arrays for hue, saturation, brightness, DGCI
hue = zeros(num_points, 1);
sat = zeros(num_points, 1);
bri = zeros(num_points, 1);
dggi = zeros(num_points, 1);

%% ========================================================================
% Step 1: Compute DGCI for each point in the point cloud
% DGCI = Dark Green Color Index
% Formula: DGCI = ((H-60)/60 + (1-S) + (1-B)) / 3
%% ========================================================================
parfor i = 1:num_points
    % Determine the dominant RGB channel
    [~, idx_max] = max(color(i,:));
    R = color(i,1); G = color(i,2); B = color(i,3);
    max_rgb = max(color(i,:));
    min_rgb = min(color(i,:));
    
    % Compute Hue
    if idx_max == 1
        hue(i) = 60 * ((G - B) / (max_rgb - min_rgb));
    elseif idx_max == 2
        hue(i) = 60 * (2 + (B - R) / (max_rgb - min_rgb));
    else
        hue(i) = 60 * (4 + (R - G) / (max_rgb - min_rgb));
    end
    
    % Compute Saturation
    sat(i) = (max_rgb - min_rgb) / max_rgb;
    
    % Compute Brightness
    bri(i) = max_rgb / 255;
    
    % Compute DGCI
    dggi(i) = ((hue(i) - 60) / 60 + (1 - sat(i)) + (1 - bri(i))) / 3;
end

%% ========================================================================
% Step 2: Compute SPAD values for each point using linear model
% SPAD = a * DGCI - b
%% ========================================================================
spad_all = a .* dggi - b;

%% ========================================================================
% Step 3: Filter noisy points
% Remove points with DGCI outside valid range (0.2 to 1)
%% ========================================================================
idx_noise = find(dggi >= 1 | dggi < 0.2);  
idx_all_error = [];

% Find faces containing noisy points
for i = 1:length(idx_noise)
    [idx_r, ~] = find(tri == idx_noise(i));
    idx_all_error = [idx_all_error; idx_r];
end
idx_all_error = unique(idx_all_error);

% Create a mask to keep only valid faces
idx_t = ones(size(tri,1),1);
idx_t(idx_all_error) = 0;

% Filter faces
plant_face = tri(logical(idx_t), :);

%% ========================================================================
% Step 4: Compute median SPAD for each valid face
% Median of SPAD values of the three vertices of each face
%% ========================================================================
leaf_spad_mesh = median([spad_all(plant_face(:,1)), ...
                         spad_all(plant_face(:,2)), ...
                         spad_all(plant_face(:,3))], 2);

end
