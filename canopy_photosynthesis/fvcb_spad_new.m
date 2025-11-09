function [A_net,A_net2] = fvcb_spad_new(leaf_spad_mesh, light_ratio_all, type, T_air)
%% FVCB_SPAD_NEW
% Compute canopy photosynthesis based on leaf SPAD values and light conditions.
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
% Inputs:
%   - leaf_spad_mesh: vector of leaf SPAD values for each leaf/facet
%   - light_ratio_all: vector of light intensity for each leaf/facet
%   - type: 'B73' or 'W64A' maize genotype
%   - T_air: air temperature (°C)
% Outputs:
%   - A_net: net CO2 assimilation for each leaf/facet under steady-state
%   - A_net2: net CO2 assimilation under non-temperature corrected condition

%% 1. Convert SPAD values to photosynthetic parameters based on genotype
if strcmp(type,'B73')
    gbs_mesh = (0.1286 .* leaf_spad_mesh - 4.8798) ./ 1000;  % bundle sheath conductance, mol/m2/s/bar
    vcmax_mesh = 0.7505 .* leaf_spad_mesh + 3.1741;           % Rubisco max rate
    vpmax_mesh = 1.7427 .* leaf_spad_mesh + 12.783;          % PEP carboxylation max rate
    jmax_mesh = 6.8604 .* leaf_spad_mesh - 119.77;           % electron transport max
    kll_mesh = -0.0002 .* leaf_spad_mesh + 0.3009;           % leaf light limitation coefficient
    Rd_mesh = -0.008 .* leaf_spad_mesh + 1.7156;             % daytime respiration
    g0_mesh = -0.0038 .* leaf_spad_mesh + 0.2374;            % stomatal intercept
    g1_mesh = 0.0819 .* leaf_spad_mesh - 0.6136;             % stomatal slope
elseif strcmp(type,'W64A')
    %% W64A genotype
    gbs_mesh = (0.0135 .* leaf_spad_mesh - 0.3849) ./ 100;   
    vcmax_mesh = 0.8881 .* leaf_spad_mesh - 0.6298;
    vpmax_mesh = 3.5927 .* leaf_spad_mesh - 20.203;
    jmax_mesh = 10.602 .* leaf_spad_mesh - 238.23;
    kll_mesh = 0.0048 .* leaf_spad_mesh + 0.089;
    Rd_mesh = 0.0408 .* leaf_spad_mesh - 0.3104;   
    g0_mesh = -0.0011 .* leaf_spad_mesh + 0.0529;
    g1_mesh = 0.0332 .* leaf_spad_mesh + 2.047;
end

%% 2. Apply minimum limits to key parameters to avoid non-physical values
jmax_min = 50;      % minimum electron transport rate
gbs_min = 0.10/1000; % minimum bundle sheath conductance
g0_min = 0.01;      % minimum stomatal intercept
gbs_mesh(gbs_mesh < gbs_min) = gbs_min;
g0_mesh(g0_mesh < g0_min) = g0_min;
jmax_mesh(jmax_mesh < jmax_min) = jmax_min;

gm = 2;  % mesophyll conductance (mol/m2/s/bar)

%% 3. Replace NaN values in light input with zero
tf = isnan(light_ratio_all);
light_ratio_all(tf) = 0;

%% 4. Set environmental conditions
path('D:\plant_model\c4_photosynthesis', path); % add model folder to MATLAB path
CO2 = 400;         % ambient CO2 (ppm)
RH = 65;           % relative humidity (%)
Oi = 0.21;         % O2 mole fraction
water_correct = 1; % water correction flag
theta_e = 0.75;    % curvature factor for light response

%% 5. Compute leaf photosynthesis for each leaf/facet using parfor
tic
parfor i = 1:length(jmax_mesh)
    % Combine photosynthetic parameters for this leaf/facet
    para = [jmax_mesh(i,1), kll_mesh(i,1), theta_e, Rd_mesh(i,1), vcmax_mesh(i,1), vpmax_mesh(i,1), gbs_mesh(i,1), gm];
    
    % Steady-state photosynthesis with temperature correction
    [A_net(i,1), gs(i,1), gv(i,1), T_leaf(i,1), E(i,1), Ci(i,1), para_temp{i,1}] = ...
        C4_steady_state_photosynthesis(light_ratio_all(i), CO2, Oi, T_air, RH, para, g0_mesh(i,1), g1_mesh(i,1), water_correct);
    
    % Non-temperature corrected photosynthesis
    [A_net2(i,1), gs2(i,1), gv2(i,1), T_leaf2(i,1), E2(i,1), Ci2(i,1), para_temp2{i,1}] = ...
        C4_steady_state_photosynthesis_NT(light_ratio_all(i), CO2, Oi, T_air, RH, para, g0_mesh(i,1), g1_mesh(i,1), water_correct);
end
toc

end
