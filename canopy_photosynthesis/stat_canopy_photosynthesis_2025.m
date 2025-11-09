function  [A_full_sum, A_neq_sum, A_aq_sum, A_aq_neq_sum, A_l_sum, A_l_neq_sum, A_detail, facet_area] = ...
    stat_canopy_photosynthesis_2025(result_str, type, T_air)
%% =========================================================================
% Function: stat_canopy_photosynthesis_2025
% Purpose : Compute canopy photosynthesis under different assumptions (full, nonsteady, acclimated)
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
% Inputs  :
%   result_str - path to the canopy simulation result file (TXT/CSV)
%   type       - plant type: 'B73' or 'W64A'
%   T_air      - air temperature (°C)
% Outputs :
%   A_full_sum      - total canopy photosynthesis (full model)
%   A_neq_sum       - total canopy photosynthesis (nonsteady model)
%   A_aq_sum        - total canopy photosynthesis using acclimated SPAD
%   A_aq_neq_sum    - total canopy photosynthesis (nonsteady, acclimated)
%   A_l_sum         - leaf-layer based photosynthesis sum
%   A_l_neq_sum     - leaf-layer nonsteady sum
%   A_detail        - detailed per-facet photosynthesis arrays
%   facet_area      - area of each canopy facet (m^2)
%% =========================================================================

%% ---------------- Step 1: Load canopy result data ----------------
result = readmatrix(result_str);           % Read simulation results
pos = result(:,6:14);                       % 3D coordinates of leaf facets
plant_idx = result(:,1);                    % Plant index for each facet
spad = result(:,15);                        % SPAD values for each facet

% Combine all points for trisurf plotting
pos_all = [pos(:,1:3); pos(:,4:6); pos(:,7:9)]; 
faces = [(1:length(pos))', (length(pos)+1:length(pos)*2)', (length(pos)*2+1:length(pos)*3)'];

facet_area = result(:,18) / 10000;         % Convert from cm^2 to m^2

% Extract PPFD (photosynthetic photon flux density) for each time step
for i = 1:15
    num = i - 1;
    PPFD(:,i) = result(:,25 + num*7);
end

leaf_spad_mesh = spad;  % Use SPAD for each leaf facet

%% ---------------- Step 2: Convert SPAD to photosynthetic parameters ----------------
if strcmp(type,'B73')
    % B73 maize genotype
    gbs_mesh   = (0.1286 .* leaf_spad_mesh - 4.8798) ./ 1000;
    vcmax_mesh = 0.7505 .* leaf_spad_mesh + 3.1741;
    vpmax_mesh = 1.7427 .* leaf_spad_mesh + 12.783;
    jmax_mesh  = 6.8604 .* leaf_spad_mesh - 119.77;
    kll_mesh   = -0.0002 .* leaf_spad_mesh + 0.3009;
    Rd_mesh    = -0.008 .* leaf_spad_mesh + 1.7156;
    g0_mesh    = -0.0038 .* leaf_spad_mesh + 0.2374;
    g1_mesh    = 0.0819 .* leaf_spad_mesh - 0.6136;
elseif strcmp(type,'W64A')
    % W64A maize genotype
    gbs_mesh   = (0.0135 .* leaf_spad_mesh - 0.3849) ./ 100;
    vcmax_mesh = 0.8881 .* leaf_spad_mesh - 0.6298;
    vpmax_mesh = 3.5927 .* leaf_spad_mesh - 20.203;
    jmax_mesh  = 10.602 .* leaf_spad_mesh - 238.23;
    kll_mesh   = 0.0048 .* leaf_spad_mesh + 0.089;
    Rd_mesh    = 0.0408 .* leaf_spad_mesh - 0.3104;
    g0_mesh    = -0.0011 .* leaf_spad_mesh + 0.0529;
    g1_mesh    = 0.0332 .* leaf_spad_mesh + 2.047;
end

% ---------------- Step 2a: Set minimum thresholds ----------------
jmax_min = 50;      % Minimum Jmax
gbs_min  = 0.10/1000;
g0_min   = 0.01;

gbs_mesh(gbs_mesh < gbs_min) = gbs_min;
g0_mesh(g0_mesh < g0_min)    = g0_min;
jmax_mesh(jmax_mesh < jmax_min) = jmax_min;

gm = 2;  % Mesophyll conductance (m/s)

%% ---------------- Step 3: Photosynthesis simulation ----------------
path('D:\plant_model\c4_photosynthesis', path);  % Add function path
CO2 = 400;        % Atmospheric CO2 (ppm)
RH  = 65;         % Relative humidity (%)
water_correct = 1;
Oi = 0.21;        % O2 concentration
theta_e = 0.75;   % Curvature parameter

% Loop over each time step
tic
for j = 1:size(PPFD,2)
    parfor i = 1:length(jmax_mesh)
        % Create parameter vector for C4 model
        para = [jmax_mesh(i), kll_mesh(i), theta_e, Rd_mesh(i), vcmax_mesh(i), vpmax_mesh(i), gbs_mesh(i), gm];
        
        % Steady-state photosynthesis
        [A_net(i,j), gs(i,j), gv(i,j), T_leaf(i,j), E(i,j), Ci(i,j), para_temp{i,j}] = ...
            C4_steady_state_photosynthesis(PPFD(i,j), CO2, Oi, T_air, RH, para, g0_mesh(i), g1_mesh(i), water_correct);
        
        % Nonsteady-state photosynthesis
        [A_net2(i,j), gs2(i,j), gv2(i,j), T_leaf2(i,j), E2(i,j), Ci2(i,j), para_temp2{i,j}] = ...
            C4_steady_state_photosynthesis_NT(PPFD(i,j), CO2, Oi, T_air, RH, para, g0_mesh(i), g1_mesh(i), water_correct);
    end
end
toc

%% ---------------- Step 4: Acclimated SPAD photosynthesis ----------------
spad_set = ones(length(leaf_spad_mesh),1) .* mean(leaf_spad_mesh);  % Mean SPAD for all leaves
for i = 1:size(PPFD,2)
    [A_mc(:,i), A_mc_NT(:,i)] = fvcb_spad_new(spad_set, PPFD(:,i), type, T_air);
end

%% ---------------- Step 5: Layer-based SPAD distribution ----------------
plant_points_t = pos_all / 100;  % Convert from cm to m
plant_face = faces;
leaf_area = LA(plant_points_t, plant_face, 1);  % Compute facet areas
total_spad = sum(leaf_area .* leaf_spad_mesh);  % Total SPAD in canopy

% Distribute SPAD across layers
[~, ~, spad_mesh] = distribute_spad_mesh_by_layer_decay(plant_points_t, plant_face, 100, 0.1, total_spad, 0.26);
spad_mesh_layer = spad_mesh ./ leaf_area;

for i = 1:size(PPFD,2)
    [A_mc_l(:,i), A_mc_l_nt(:,i)] = fvcb_spad_new(spad_mesh_layer, PPFD(:,i), type, T_air);
end

%% ---------------- Step 6: Compute per-facet photosynthesis ----------------
for i = 1:size(PPFD,2)
    leaf_A_full(:,i)      = A_net(:,i) .* facet_area;
    leaf_A_neq(:,i)       = A_net2(:,i) .* facet_area;
    leaf_A_aq(:,i)        = A_mc(:,i) .* facet_area;
    leaf_A_aq_neq(:,i)    = A_mc_NT(:,i) .* facet_area;
    leaf_A_l(:,i)         = A_mc_l(:,i) .* facet_area;
    leaf_A_l_neq(:,i)     = A_mc_l_nt(:,i) .* facet_area;
end

%% ---------------- Step 7: Sum over canopy ----------------
A_full_sum    = sum(leaf_A_full,1);
A_neq_sum     = sum(leaf_A_neq,1,'omitnan');
A_aq_sum      = sum(leaf_A_aq,1);
A_aq_neq_sum  = sum(leaf_A_aq_neq,1);
A_l_sum       = sum(leaf_A_l,1,'omitnan');
A_l_neq_sum   = sum(leaf_A_l_neq,1,'omitnan');

A_detail = {leaf_A_full, leaf_A_neq, leaf_A_aq, leaf_A_aq_neq, leaf_A_l, leaf_A_l_neq};

end
