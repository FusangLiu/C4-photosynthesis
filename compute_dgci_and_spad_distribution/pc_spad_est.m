%% ========================================================================
%  Script Title : Compute DGCI and Estimate SPAD from Plant Point Cloud
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
%  Date         : November 2023
%
%  Description :
%  This script reads a 3D reconstructed plant model (from Geomagic output),
%  computes the Dark Green Color Index (DGCI) for each point in the cloud
%  based on RGB color values, and estimates SPAD (chlorophyll content)
%  using a linear relationship with DGCI.
%
%  It provides 3D visualizations of:
%     1. DGCI distribution across the plant canopy
%     2. SPAD distribution (chlorophyll level)
%
%  Requirements :
%     -  .ply files exported from Geomagic:
%         1) "nerf_b73-pc2-meshed_plant.ply"  → triangle mesh (surface)
%         2) "nerf_b73-pc2-meshed_plant.ply" → colorized point cloud
%
%  Output :
%     - 3D scatter plots of DGCI and SPAD
%     - (Optional) high-resolution SPAD image can be exported as PNG
% ========================================================================


%% ========================================================================
%  Step 1: Read the plant model from Geomagic files
% ========================================================================
[tri, ~] = plyread('nerf_b73-pc2-meshed_plant.ply', 'tri');  % Triangular mesh data
pc = pcread('nerf_b73-pc2-meshed_plant.ply');             % Point cloud data with RGB
pos = pc.Location;                                            % Nx3 coordinates

% Display semi-transparent plant mesh
figure;
trisurf(tri, pos(:,1), pos(:,2), pos(:,3), ...
    'FaceColor', [0.133, 0.545, 0.133], ...  % RGB: green
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'none');
axis off;
axis equal;
set(gcf, 'Color', [1,1,1]);


%% ========================================================================
%  Step 2: Compute DGCI (Dark Green Color Index)
% ========================================================================
color = double(pc.Color);     % RGB color values (Nx3, range 0–255)
nPoints = size(color, 1);

% Preallocate arrays for speed
hue = zeros(nPoints, 1);
sat = zeros(nPoints, 1);
bri = zeros(nPoints, 1);
dgci = zeros(nPoints, 1);

% Compute HUE, SATURATION, BRIGHTNESS, and DGCI in parallel
parfor i = 1:nPoints
    pic_R = color(i,1);
    pic_G = color(i,2);
    pic_B = color(i,3);

    max_rgb = max([pic_R, pic_G, pic_B]);
    min_rgb = min([pic_R, pic_G, pic_B]);

    % Avoid division by zero
    if max_rgb == min_rgb
        hue(i,1) = 0;
    else
        % Determine hue angle based on RGB dominance
        [~, idx] = max(color(i,:));
        if idx == 1
            hue(i,1) = 60 * ((pic_G - pic_B) / (max_rgb - min_rgb));
        elseif idx == 2
            hue(i,1) = 60 * (2 + (pic_B - pic_R) / (max_rgb - min_rgb));
        else
            hue(i,1) = 60 * (4 + (pic_R - pic_G) / (max_rgb - min_rgb));
        end
    end

    % Compute saturation and brightness
    sat(i,1) = (max_rgb - min_rgb) / max_rgb;
    bri(i,1) = max_rgb / 255;

    % Dark Green Color Index (DGCI)
    % (This version follows the normalized composite form)
    dgci(i,1) = ((hue(i,1) - 60) / 60 + (1 - sat(i,1)) + (1 - bri(i,1))) / 3;
end


%% ========================================================================
%  Step 3: Visualize DGCI distribution
% ========================================================================
valid_idx = find(dgci < 1);   % Filter valid range

figure;
scatter3(pos(valid_idx,1), pos(valid_idx,2), pos(valid_idx,3), ...
    ones(length(valid_idx),1), dgci(valid_idx,1), 'filled');
axis equal; axis off;
colorbar;
colormap('turbo');
caxis([0 1]);
set(gcf, 'Color', [1 1 1]);
title('DGCI Distribution (Dark Green Color Index)', 'FontSize', 12);
view(82, 27);


%% ========================================================================
%  Step 4: Estimate and visualize SPAD (chlorophyll index)
% ========================================================================
% Linear model derived from calibration: SPAD = a * DGCI + b
spad = 84.47 .* dgci(valid_idx,1) - 2.487;

figure;
scatter3(pos(valid_idx,1), pos(valid_idx,2), pos(valid_idx,3), ...
    ones(length(valid_idx),1), spad, 'filled');
axis equal; axis off;
colorbar;
caxis([0 prctile(spad, 99.9)]);  % reduce outlier effect
colormap(green_color_bar);        % custom colormap (assumed predefined)
set(gcf, 'Color', [1 1 1]);
view(82, 27);
title('SPAD Distribution Estimated from DGCI', 'FontSize', 12);

% Optionally save figure (uncomment to export)
% print(gcf, '-dpng', '-r400', 'spad.png');
