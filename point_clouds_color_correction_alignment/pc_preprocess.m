%% ========================================================================
%  Script Title : 3D Point Cloud Color Correction and Alignment Pipeline
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
%  Date         : November 2023
%
%  Description  :
%  This MATLAB script performs a full pipeline for color correction,
%  noise filtering, segmentation, plane fitting, alignment, and scaling
%  of plant 3D point clouds captured with RGB sensors.
%
%  Workflow Steps :
%   1. Load and visualize raw point cloud
%   2. Perform color correction using a color checker chart
%   3. Apply color correction to the 3D point cloud
%   4. Filter noise and separate background/foreground
%   5. Fit a plane for geometric alignment
%   6. Align plant to reference plane
%   7. Apply scale transformation
%   8. Save the corrected point cloud
%
%  Note:
%  - Steps 2 and 7 require manual user input (selecting points in figure).
%  - All intermediate visualizations are retained for verification.
% ========================================================================


%% ===============================================================
%  Step 1: Load and visualize original point cloud
% ================================================================
clear; close all;
str = 'nerf_b73-pc2';  % base filename (modify as needed)

% Read and display the raw point cloud
pc_d = pcread([str, '.ply']);
figure;
pcshow(pc_d);
axis off; axis equal;

% Save visualization as PNG for color correction
print(gcf, '-dpng', '-r400', [str, '.png']);


%% ===============================================================
%  Step 2: Perform color correction using color checker
%  NOTE: This step requires manual selection of reference points.
% ================================================================

I = imread([str, '.png']);
figure; imshow(I);

% --- MANUAL OPERATION REQUIRED ---
% Select the following points on the image:
%  1. black point
%  2. white point
%  3. dark skin tone point
%  4. bluish-green point
blackPoint = drawpoint;
whitePoint = drawpoint;
darkSkinPoint = drawpoint;
bluishGreenPoint = drawpoint;

cornerPoints = [blackPoint.Position;
                whitePoint.Position;
                darkSkinPoint.Position;
                bluishGreenPoint.Position];

% Detect color checker chart
chart = colorChecker(I, "RegistrationPoints", cornerPoints);
figure; displayChart(chart);

% Measure color and compute color correction matrix (CCM)
[colorTable, ccm] = measureColor(chart);
figure; displayColorPatch(colorTable);

% Apply CCM to the original image
I_cc = imapplymatrix(ccm(1:3,:)', I, ccm(4,:));
figure; imshow(I_cc);
title('Color-Corrected Image (1st Pass)');

% Recalculate correction using the color-corrected image
chart_cc = colorChecker(I_cc, "RegistrationPoints", cornerPoints);
colorTable_cc = measureColor(chart_cc);
figure; displayColorPatch(colorTable_cc);

% Second-pass correction for improved accuracy
[colorTable2, ccm2] = measureColor(chart_cc);
figure; displayColorPatch(colorTable2);
I_cc2 = imapplymatrix(ccm2(1:3,:)', I_cc, ccm2(4,:));
figure; imshow(I_cc2);
title('Color-Corrected Image (2nd Pass)');

% Display final color patch chart
chart_cc2 = colorChecker(I_cc2, "RegistrationPoints", cornerPoints);
colorTable_cc2 = measureColor(chart_cc2);
figure; displayColorPatch(colorTable_cc2);


%% ===============================================================
%  Step 3: Apply color correction matrix to 3D point cloud
% ================================================================

pos = pc_d.Location;
color_o = double(pc_d.Color);

% Apply color correction to RGB values
color = [color_o, ones(length(color_o),1)] * ccm;
pc_t = pointCloud(pos, 'Color', uint8(color));

% Display and save the color-corrected point cloud
figure; pcshow(pc_t);
pcwrite(pc_t, [str, '_cc.ply']);


%% ===============================================================
%  Step 4: Noise filtering and background/foreground separation
%  NOTE: Requires background.ply and foreground.ply as references.
% ================================================================

% Denoise the corrected point cloud
pc = pcdenoise(pcread([str, '_cc.ply']));
figure; pcshow(pc);

% Load background and foreground reference clouds
pc_b = pcread('background.ply');
pc_f = pcread('foreground.ply');

% Perform SVM-based color segmentation (custom function)
[~, pc_d, ~, background] = SVM_color_pointcloud_seg2(pc, pc_b, pc_f, 0, [], 1);

veg = pc_d;  % keep only vegetation points


%% ===============================================================
%  Step 5: Plane fitting for background and orientation alignment
% ================================================================

panel_pc = background;

% Fit a plane to the background (e.g., calibration panel or table)
[model, inlierIndices, ~] = pcfitplane(panel_pc, 0.1);
plane1 = select(panel_pc, inlierIndices);
plane_center = mean(plane1.Location);
Parameters = model.Parameters;

% Display fitted plane and its normal vector
figure;
pcshow(plane1);
hold on;
quiver3(plane_center(1), plane_center(2), plane_center(3), ...
        Parameters(1), Parameters(2), Parameters(3));


%% ===============================================================
%  Step 6: Align plant point cloud relative to the fitted plane
% ================================================================

plant_pc = veg;
pos = plant_pc.Location;
pos_plane = plane1.Location;

% Define normal vector of fitted plane
normal_pnts = Parameters(:, 1:3)';

% Compute rotation to align with Z-axis
n0 = [0, 0, 1];
rotation_axis = cross(normal_pnts', n0);
theta = acos((n0 * normal_pnts) / (norm(n0) * norm(normal_pnts')));
Rv = rotation_axis / norm(rotation_axis) * theta;
rotationMatrix = rotationVectorToMatrix(Rv);

% Extend to 4x4 homogeneous matrix
rotationMatrix(:,4) = 0;
rotationMatrix(4,:) = [0 0 0 1];

% Apply rotation to plant and plane points
pos_t = [pos, ones(size(pos,1),1)];
pos_tp = [pos_plane, ones(size(pos_plane,1),1)];
pos_end = pos_t * rotationMatrix;
pos_tp_end = pos_tp * rotationMatrix;

% Center the rotated point cloud
pos_end = pos_end(:,1:3) - mean(pos_tp_end(:,1:3));
pos_tr = pointCloud(pos_end, 'Color', plant_pc.Color);
figure; pcshow(pos_tr);


%% ===============================================================
%  Step 7: Scale the aligned point cloud
%  NOTE: Manual cursor measurement required.
% ================================================================

% --- MANUAL OPERATION REQUIRED ---
% Define cursor_info and cursor_info1 by clicking two points 
% whose physical distance is known (e.g., 12 mm).
scale = 0.012 / sqrt(sum((cursor_info.Position - cursor_info1.Position).^2));

% Apply scaling
pos_all = pos_end * scale;
pos_tr = pointCloud(pos_all, 'Color', plant_pc.Color);


%% ===============================================================
%  Step 8: Save and display final corrected 3D model
% ================================================================

figure;
pcshow(pos_tr);
xlabel('X'); ylabel('Y'); zlabel('Z');
view(0,0); axis equal;

% Save final corrected point cloud
pcwrite(pos_tr, [str, '_cc_o.ply']);
