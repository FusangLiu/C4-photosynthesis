function [pL, pCloud, model, pCloud2] = SVM_color_pointcloud_seg2(p, background, foreground, method, models, show)
%% ========================================================================
%  Function Title : SVM-based Point Cloud Color Segmentation
%  Author         : Fusang Liu
%  Affiliation    : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                   Wageningen University & Research (WUR)
%  Date           : November 2021
%
%  Description :
%  This function segments a 3D point cloud (p) into foreground and
%  background regions based on RGB color information using either an
%  SVM (Support Vector Machine) or Random Forest (TreeBagger) classifier.
%
%  INPUTS :
%   p           - input point cloud object to be segmented
%   background  - point cloud of known background color samples
%   foreground  - point cloud of known foreground (vegetation) samples
%   method      - classifier selection:
%                   0 = SVM (default)
%                   1 = TreeBagger (Random Forest)
%   models      - pre-trained model (optional; empty if training new)
%   show        - logical flag (1 = show results, 0 = silent)
%
%  OUTPUTS :
%   pL          - Nx3 coordinates of segmented (foreground) points
%   pCloud      - pointCloud object of segmented foreground
%   model       - trained classifier (SVM or TreeBagger)
%   pCloud2     - pointCloud object of removed background
%
%  Usage Example :
%   [pL, pCloud, model, pCloud2] = SVM_color_pointcloud_seg2(pc, pc_b, pc_f, 0, [], 1);
%
%  Dependencies :
%   - Requires MATLAB’s Statistics and Machine Learning Toolbox
%   - Designed to be used with RGB point clouds
% ========================================================================

%% ===============================================================
%  Step 1: Extract RGB color data from the input point cloud
% ================================================================
pic = p.Color;   % RGB color matrix (Nx3)

%% ===============================================================
%  Step 2: Train a classification model (if not provided)
% ================================================================
if isempty(models)
    % Extract color samples from training point clouds
    TrainData_background = background.Color;       
    TrainData_foreground = foreground.Color;

    % Assign class labels: 0 = background, 1 = foreground
    TrainLabel = [zeros(length(TrainData_background),1); ...
                  ones(length(TrainData_foreground),1)];

    % Combine training data and convert to double precision
    TrainData = double([TrainData_background; TrainData_foreground]);

    % Choose classifier type
    if method == 1
        % Random Forest classifier
        model = TreeBagger(600, TrainData, TrainLabel);
    else
        % Support Vector Machine classifier
        model = fitcsvm(TrainData, TrainLabel);
    end
else
    % Use pre-trained model
    model = models;
end


%% ===============================================================
%  Step 3: Predict class labels for all points in the input cloud
% ================================================================
TestData = double(pic);
TestLabel = predict(model, TestData);   % Predicted labels (0 or 1)

% Convert single label vector to logical RGB mask (Nx3)
ind = [TestLabel, TestLabel, TestLabel];
ind = logical(ind);

% Create a color mask for foreground points
pic_seg = pic;
pic_seg(~ind) = 0;  % Set background points to black (0,0,0)
p.Color = pic_seg;  % Update point cloud colors


%% ===============================================================
%  Step 4: Separate foreground and background points
% ================================================================
ps = double(pic_seg);

% Combine positions and colors
a = [p.Location, ps];  % all points
% Keep only non-black points (foreground)
a(any(a(:,4:6)==0,2),:) = [];

% Create another matrix for background
b = [p.Location, ps];
% Keep only black (background) points
b(any(b(:,4:6)~=0,2),:) = [];

% Extract position and color data for both classes
pL  = single(a(:,1:3));    % foreground coordinates
pCr = uint8(a(:,4:6));     % foreground colors
pL2 = single(b(:,1:3));    % background coordinates
pCr2 = uint8(b(:,4:6));    % background colors

% Build separate point cloud objects
pCloud  = pointCloud(pL,  'Color', pCr);   % Foreground
pCloud2 = pointCloud(pL2, 'Color', pCr2);  % Background


%% ===============================================================
%  Step 5: Optional visualization of segmentation result
% ================================================================
if show
    % Display segmented foreground
    figure;
    pcshow(pCloud);
    axis equal; axis off;
    set(gcf, 'Color', [1,1,1]);

    % Display background points
    figure;
    pcshow(pCloud2);
    axis equal; axis off;
    set(gcf, 'Color', [1,1,1]);
end

end