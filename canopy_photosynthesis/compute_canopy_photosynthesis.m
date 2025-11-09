%% ========================================================================
%  Script Title : Canopy Photosynthesis Simulation for B73 and W64A
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
%  Date         : November 2024
%
%  Description :
%  This script:
%     1. Reads 3D plant point clouds and meshes from Geomagic output
%     2. Computes SPAD values for each plant facet
%     3. Constructs canopy models with randomized planting positions
%     4. Performs ray tracing to simulate light interception
%     5. Aggregates canopy photosynthesis results
%
%  Requirements :
%     - Plant mesh and point cloud .ply files
%     - Random matrices: 'random_mat_4_4.mat'
%     - Custom functions: compute_spad_pc, canopy_model_spad, stat_canopy_photosynthesis_2025, fvcb_spad_new
% ========================================================================

clear; close all;

%% ========================================================================
% Step 1: Set SPAD calibration parameters for B73 jointing stage
% ========================================================================
a = 160.69; b = 26.946;    % linear SPAD model: SPAD = a * DGCI + b
% a=85.734; b=3.8895; %w64a nerf seeding
% a=105.69; b=15.965;%w64a nerf jointing
% a=184.78; b=45.195; %b73 nerf seeding
out_name = 'D:\plant_model\c4_photosynthesis\github\canopy_photosynthesis\group_model\hd15_b73_jointing_r';

%% ========================================================================
% Step 2: Read and process plant point clouds in center of canopy
% ========================================================================
plant_pts = cell(6,1);
plant_spad = cell(6,1);
plant_faces = cell(6,1);
plant_idx = cell(6,1);

for i = 1:3
    num = i;
    str = ['jointing-nerf_b73_', num2str(num)];
    
    % Load plant mesh and point cloud
    [tri, ~] = plyread([str, '-meshed_plant.ply'], 'tri');
    pc = pcread([str, '-meshed_plant.ply']);
    pos = pc.Location;
    
    % Visualize plant mesh
    figure;
    trisurf(tri, pos(:,1), pos(:,2), pos(:,3), ...
        'FaceColor', [0.133, 0.545, 0.133], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    axis equal; axis off;
    set(gcf, 'Color', [1,1,1]);
    view(112.4, 9.5);
    
    % Store plant data
    plant_pts{num,1} = pos;
    [plant_spad{num,1}, plant_faces{num,1}] = compute_spad_pc(pc, tri, a, b);
    plant_idx{num,1} = ones(length(plant_faces{num,1}),1) * num;
end

%% ========================================================================
% Step 3: Read and process plant point clouds in outside of canopy
% ========================================================================
for i = 1:3
    num = i + 3;
    str = ['jointing-nerf_b73_', num2str(i)];
    
    [tri, ~] = plyread([str, '-meshed_plant.ply'], 'tri');
    pc = pcread([str, '-meshed_plant.ply']);
    pos = pc.Location;
    
    figure;
    trisurf(tri, pos(:,1), pos(:,2), pos(:,3), ...
        'FaceColor', [0.133, 0.545, 0.133], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    axis equal; axis off;
    set(gcf, 'Color', [1,1,1]);
    view(112.4, 9.5);
    
    plant_pts{num,1} = pos;
    [plant_spad{num,1}, plant_faces{num,1}] = compute_spad_pc(pc, tri, a, b);
    plant_idx{num,1} = ones(length(plant_faces{num,1}),1) * num;
end

%% ========================================================================
% Step 4: Generate random matrices for canopy modeling
% ========================================================================
load('random_mat_4_4.mat'); % load pre-defined random matrix

all_rm = cell(3,1);
all_dir = cell(3,1);

% Generate random positioning and direction matrices
for i = 1:3
    gourp_rm = randi([4,6], 9,7);
    center_rm = randi([1,3], 3,3);
    gourp_rm(4:6,3:5) = center_rm;
    all_rm{i,1} = gourp_rm;
    
    gourp_dir = rand(9,7);
    all_dir{i,1} = gourp_dir;
end

%% ========================================================================
% Step 5: Build canopy model and compute SPAD for each facet
% ========================================================================
plant_dis = 0.16; % distance between plants (m)
row_dis = 0.41;   % distance between rows (m)
plant_num = 4;    % plants per row
row_num = 4;      % number of rows

group_rp_all = cell(3,1);
group_re_all = cell(3,1);
group_spad_all = cell(3,1);
canopy_idx_all_g = cell(3,1);

for i = 1:3
    [group_rp_all{i,1}, group_re_all{i,1}, group_spad_all{i,1}, canopy_idx_all_g{i,1}] = ...
        canopy_model_spad(plant_faces, plant_pts, plant_spad, plant_idx, ...
        plant_num, row_num, plant_dis, row_dis, all_rm{i,1}, all_dir{i,1});
end

%% ========================================================================
% Step 6: Visualize canopy SPAD and save data
% ========================================================================
load('green_color_bar.mat'); % custom green colormap

for i = 1:3
    group_rp = group_rp_all{i,1};
    group_re = group_re_all{i,1};
    group_spad = group_spad_all{i,1};
    canopy_idx_all = canopy_idx_all_g{i,1};
    
    figure;
    trisurf(group_re, group_rp(:,1), group_rp(:,2), group_rp(:,3), group_spad, ...
        'FaceAlpha', 0.9, 'EdgeColor', 'none');
    caxis([prctile(group_spad,1) prctile(group_spad,95)]);
    axis equal
    set(gcf,'Color',[1,1,1]);
    colormap(green_color_bar);
    colorbar;
    view(82,27);
    
    % Save facet information to text file for ray tracing
    canopy_facet_num = length(canopy_idx_all);
    input_mat = [canopy_idx_all(:,1), ones(canopy_facet_num,1), ...
                 canopy_idx_all(:,2), zeros(canopy_facet_num,1), ...
                 zeros(canopy_facet_num,1), ...
                 group_rp(group_re(:,1),:).*100, group_rp(group_re(:,2),:).*100, group_rp(group_re(:,3),:).*100, ...
                 group_spad, ones(canopy_facet_num,1).*0.075, ones(canopy_facet_num,1).*0.075];
    writematrix(input_mat, [out_name,num2str(i),'_new.txt'], 'Delimiter',' ');
end

%% ========================================================================
% Step 7: Ray tracing simulation using fastTracer
% ========================================================================
path_rt = 'D:\plant_model\c4_photosynthesis\github\canopy_photosynthesis';
path_in = 'D:\plant_model\c4_photosynthesis\github\canopy_photosynthesis\group_model\';
path_out = 'D:\plant_model\c4_photosynthesis\github\canopy_photosynthesis\canopy_photosyntheis_result\';

files = dir([path_in, '*txt']);
parfor i = 1:3
    file_name = files(i).name;
    out_name_rt = split(file_name, ".");
    cmd = [path_rt,'\fastTracerV1.22.exe -D 24 56 61.5 143.5 0 230 -L 31 -S 12 -A 0.7 -d 120 -W 5 1 19 -n 0.1 -m ', ...
            path_in, file_name, ' -o ', path_out, out_name_rt{1}, '_result.txt -z 0.5'];
    system(cmd);
end

%% ========================================================================
% Step 8: Read ray tracing results and compute canopy photosynthesis
% ========================================================================
files_result = dir([path_out, '*txt']);

for i = 1:3
    [A_full_sum(i,:), A_neq_sum(i,:), A_aq_sum(i,:), A_l_sum(i,:), A_detail{i}, facet_area{i}] = ...
        stat_canopy_photosynthesis_2025([path_out, files_result(i).name], 'B73',28);
end


%% ========================================================================
% Step 9: Post-process and reshape hourly and daily results
% ========================================================================
% Example of hourly energy per leaf, summed over canopy
for i=1:3
    A_full_sum_n(i,:)=[min(A_full_sum(i,:)).*ones(1,4),A_full_sum(i,:),min(A_full_sum(i,:)).*ones(1,5)];
    A_neq_sum_n(i,:)=[min(A_neq_sum(i,:)).*ones(1,4),A_neq_sum(i,:),min(A_neq_sum(i,:)).*ones(1,5)];
    A_aq_sum_n(i,:)=[min(A_aq_sum(i,:)).*ones(1,4),A_aq_sum(i,:),min(A_aq_sum(i,:)).*ones(1,5)];
end

afull=reshape(A_full_sum_n',3*24,1)/(0.4*1);
aneq=reshape(A_neq_sum_n',3*24,1)/(0.4*1);
aq=reshape(A_aq_sum_n',3*24,1)/(0.4*1);

time=[1:24]';
time_all=repmat(time,3,1);
all_result=[afull,aneq,aq,time_all];


for i=1:3
all_day_afull(i,1)=hour_to_day_energy(A_full_sum_n(i,:),24);
all_day_aneq(i,1)=hour_to_day_energy(A_neq_sum_n(i,:),24);
all_day_aq(i,1)=hour_to_day_energy(A_aq_sum_n(i,:),24);
end

all_day_result=[all_day_afull,all_day_aneq,all_day_aq]/10^6;
