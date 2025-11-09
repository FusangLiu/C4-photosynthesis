% ============================================================
% Description: C4 photosynthesis model parameter estimation
% Author:      Fusang Liu
% Date:        2022-04-26
% ============================================================

close all
clear

% ------------------------------------------------------------
% Add paths for MCMC diagnostic tools and core functions
% ------------------------------------------------------------
addpath('D:\plant_model\c4_photosynthesis\matlab-bayesian-estimation-master\fileExchange\mcmcdiag');
addpath('D:\plant_model\c4_photosynthesis\matlab-bayesian-estimation-master');

% ------------------------------------------------------------
% Load raw gas-exchange and fluorescence data
% ------------------------------------------------------------
name = 'D:\plant_model\c4_photosynthesis\raw_gas_exchange\raw_data_230210_2.xlsx';
A_Q_am  = xlsread(name, 'A_Q_am');
A_Ci_am = xlsread(name, 'A_Ci_am');

% A–Q data
PPFD = A_Q_am(:,1);      % Incident light (μmol m⁻² s⁻¹)
A    = A_Q_am(:,2);      % Net assimilation rate (μmol m⁻² s⁻¹)
Ci   = A_Q_am(:,3);      % Intercellular CO₂ (μmol mol⁻¹)
Y2   = A_Q_am(:,4);      % PSII quantum yield (unitless)

% A–Ci data
PPFD_ci = A_Ci_am(:,1);
A_ci    = A_Ci_am(:,2);
Ci_ci   = A_Ci_am(:,3);
Y2_ci   = A_Ci_am(:,4);

% ------------------------------------------------------------
% Quick data visualization
% ------------------------------------------------------------
figure; plot(PPFD, A, '.');     xlabel('PPFD'); ylabel('A');
figure; plot(PPFD, Y2, '.');    xlabel('PPFD'); ylabel('Y2');
figure; plot(PPFD, Ci, '.');    xlabel('PPFD'); ylabel('Ci');
figure; plot(Ci_ci, A_ci, '.'); xlabel('Ci');   ylabel('A');
figure; plot(Ci_ci, Y2_ci, '.');xlabel('Ci');   ylabel('Y2');

% ------------------------------------------------------------
% Estimate respiration in the light (R_light)
% and low-light PSII yield (Y2_LL) using Yin et al. (2011)
% ------------------------------------------------------------
low_light_idx = find(PPFD <= 150);

para_reg     = polyfit(PPFD(low_light_idx).*Y2(low_light_idx)./3, A(low_light_idx), 1);
para_reg_y2  = polyfit(PPFD(low_light_idx), Y2(low_light_idx), 1);
R_light      = -para_reg(2);
Y2_LL        = para_reg_y2(2);
s            = para_reg(1);

% Plot Y2 vs. PPFD under low light
x1 = 1:150;
y1 = polyval(para_reg_y2, x1);
figure
plot(PPFD(low_light_idx), Y2(low_light_idx), 'o')
hold on; plot(x1, y1)
xlabel('PPFD'); ylabel('Y2')
ylim([0.5, 0.75])
set(gcf,'Color',[1 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5 2.2])

% ------------------------------------------------------------
% Calculate ATP production rate (JATP) using Yin’s method
% ------------------------------------------------------------
idx_ci_am = find(Ci_ci >= 75);
PPFD_fit_all_am = [PPFD; PPFD_ci(idx_ci_am)];
Y2_fit_all_am   = [Y2;   Y2_ci(idx_ci_am)];

J_ATP_yin = (s .* Y2_fit_all_am .* PPFD_fit_all_am) ./ (1 - 0.4);

% ------------------------------------------------------------
% Fit light-response curve of ATP production (non-rectangular hyperbola)
% ------------------------------------------------------------
para0 = [s*Y2_LL/(1-0.4) 140 0.75]; % initial guess [J_ATP_LL, J_ATP_sat, theta]
[para_fit_atp, ~] = lsqcurvefit(@non_rectangular_hyperbola_curve, ...
    para0, PPFD_fit_all_am, J_ATP_yin, ...
    [s*Y2_LL/(1-0.4) 0 0.75], [s*Y2_LL/(1-0.4)+inf inf 0.75]);

J_ATP_LL_fit = para_fit_atp(1);
J_ATP_sat    = para_fit_atp(2);
theta        = para_fit_atp(3);

% Plot JATP fit
PPFD_m = 0:1:2500;
for i = 1:length(PPFD_m)
    J_ATP_m_am(i) = non_rectangular_hyperbola_curve(para_fit_atp, PPFD_m(i));
end
figure
plot(PPFD_fit_all_am, J_ATP_yin, 'o')
hold on; plot(PPFD_m, J_ATP_m_am, '-')
xlabel('PPFD'); ylabel('JATP')
set(gcf,'Color',[1 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5 2.2])

% ------------------------------------------------------------
% Estimate bundle-sheath conductance (gBS) using J/J method
% ------------------------------------------------------------
JJ_idx_AQ  = find(PPFD <= 500);
JJ_idx_Aci = find(Ci_ci >= 100);

PPFD_JJ = [PPFD(JJ_idx_AQ); PPFD_ci(JJ_idx_Aci)];
Y2_JJ   = [Y2(JJ_idx_AQ);   Y2_ci(JJ_idx_Aci)];
A_JJ    = [A(JJ_idx_AQ);    A_ci(JJ_idx_Aci)];
Ci_JJ   = [Ci(JJ_idx_AQ);   Ci_ci(JJ_idx_Aci)];

J_ATP_yin_JJ = (s .* Y2_JJ .* PPFD_JJ) ./ (1 - 0.4);

parameters_est_gbs = {A_JJ, Ci_JJ, R_light};
[gbs_jj, ~] = lsqcurvefit(@gBS_est, 1, parameters_est_gbs, J_ATP_yin_JJ, 0.1, 20);

J_ATP_bs = gBS_est(gbs_jj, parameters_est_gbs);

% Plot gBS estimation
figure
plot(PPFD_JJ, J_ATP_yin_JJ, '*')
hold on; plot(PPFD_JJ, J_ATP_bs, '.')
xlabel('PPFD'); ylabel('JATP (Bundle Sheath)')
set(gcf,'Color',[1 1 1])

% ------------------------------------------------------------
% Estimate Vpmax from low-Ci data
% ------------------------------------------------------------
vp_idx_Aci = find(Ci_ci <= 75);
A_vp  = A_ci(vp_idx_Aci);
Ci_vp = Ci_ci(vp_idx_Aci);

parameters_est_vp = {A_vp, Ci_vp, R_light, gbs_jj/1000};
[Vpmax, ~] = lsqcurvefit(@vpmax_est, 100, parameters_est_vp, A_vp, 0, 1000);

A_mod = vpmax_est(Vpmax, parameters_est_vp);
figure
plot(A_mod, A_vp, '*')
hold on; plot(A_vp, A_vp, '_')
xlabel('Modelled A'); ylabel('Measured A')

% ------------------------------------------------------------
% MCMC setup for Bayesian parameter estimation
% ------------------------------------------------------------
nSamples     = 40000; % samples per chain
nChains      = 4;     % number of MCMC chains
thinSteps    = 1;
burnInSteps  = 20000; % burn-in samples
parameters   = {'VCMAX','VPMAX','GBS'};

% ------------------------------------------------------------
% Prepare full dataset for MCMC fitting
% ------------------------------------------------------------
A_all = [A; A_ci];
data_in = [PPFD, Ci, Y2; PPFD_ci, Ci_ci, Y2_ci];
data_in(:,4) = 21 * 10000;  % constant O2 concentration

PPFD_all = data_in(:,1);
Ci_all   = data_in(:,2);
Y2_all   = data_in(:,3);
O2_all   = data_in(:,4);

Jatp_all_yin = (s .* Y2_all .* PPFD_all) ./ (1 - 0.4);
R_light_all  = ones(length(Jatp_all_yin),1) .* R_light;

% Compute non-rectangular hyperbola-based JATP
for i = 1:size(data_in,1)
    Jatp_all(i) = non_rectangular_hyperbola_curve(para_fit_atp, data_in(i,1));
end

% Create data structure for JAGS input
for i = 1:nChains
    dataList(i) = struct('Aobs',A_all,'JATP',Jatp_all_yin,'CI',Ci_all,'OM',O2_all,...
        'RD',R_light_all,'KP',80,'KMC',650,'KMO',450000,'SCO',2590.7,...
        'AL',0.15,'X',0.4,'GM',2,'numPts',size(A_all,1));
end

% Initial parameter guesses
for i = 1:nChains
    initsList(i) = struct('VCMAX',40,'VPMAX',80,'GBS',3/1000);
end

% ------------------------------------------------------------
% Run MCMC using JAGS model
% ------------------------------------------------------------
model = fullfile(pwd, 'C4_t2.txt');
[samples, stats, mcmcChain] = matjags(...
    dataList, model, initsList, ...
    'monitorparams', parameters, ...
    'nChains', nChains, ...
    'nSamples', nSamples, ...
    'nBurnin', burnInSteps, ...
    'thin', thinSteps, ...
    'doParallel', 1);

% ------------------------------------------------------------
% Process and visualize MCMC results
% ------------------------------------------------------------
mcmcChain2 = mbe_restructChains(mcmcChain);
mbe_diagMCMC(mcmcChain2);
mbe_plotPairs(mcmcChain2, 1000)

% Extract posterior means
for i=1:4
    [count, vcmax_m] = hist(mcmcChain2.VCMAX1(:,i),100); [~,idx]=max(count); vcmax_temp(i)=vcmax_m(idx);
    [count, vpmax_m] = hist(mcmcChain2.VPMAX1(:,i),100); [~,idx]=max(count); vpmax_temp(i)=vpmax_m(idx);
    [count, gbs_m]   = hist(mcmcChain2.GBS1(:,i),100);   [~,idx]=max(count); gbs_temp(i)=gbs_m(idx);
end
vcmax = mean(vcmax_temp);
vpmax = mean(vpmax_temp);
gbs   = mean(gbs_temp);
gm    = 2;

% ------------------------------------------------------------
% Compute model predictions
% ------------------------------------------------------------
for i = 1:length(Jatp_all)
    A_min2(i,1) = C4_mcmc(Jatp_all(i), Ci_all(i), O2_all(i), R_light_all(i), gbs, vcmax, vpmax, gm);
    A_min(i,1)  = C4_mcmc(Jatp_all_yin(i), Ci_all(i), O2_all(i), R_light_all(i), gbs, vcmax, vpmax, gm);
end

% Plot measured vs. modeled A
A_jatp = (1-0.4).*Jatp_all_yin/3 - R_light;
figure; plot(A_min, A_all,'.'); hold on
plot(A_min2, A_all,'.'); plot(A_all, A_all,'-')
xlabel('Calculated A'); ylabel('Measured A')

% Light response plots
figure
plot(PPFD, A,'.'); hold on
plot(PPFD, A_min(1:length(A)),'-');
plot(PPFD, A_min2(1:length(A)),'.');
xlabel('PPFD'); ylabel('A')

% CO2 response plots
data_aci = sortrows([Ci_ci, A_min(length(A)+1:end)], 1);
figure
plot(Ci_ci, A_ci,'.'); hold on
plot(data_aci(:,1), data_aci(:,2),'-');
plot(Ci_ci, A_min2(length(A)+1:end),'.');
xlabel('Ci'); ylabel('A')

% ------------------------------------------------------------
% Fit stomatal conductance (gs) using Ball-Berry model
% ------------------------------------------------------------
para0 = [0.05, 4]; 
gs    = A_Q_am(:,5);
input = A_Q_am(:,[2,6,7,8]);
[para_fit_gs, ~] = lsqcurvefit(@stomatal_conductance_BB_est, para0, input, gs, [0 0], [1 10]);

gs_fit = stomatal_conductance_BB_est(para_fit_gs, input);
figure
plot(input(:,1), gs_fit,'o'); hold on; plot(input(:,1), gs,'*')
xlabel('A'); ylabel('gs')

% ------------------------------------------------------------
% Fit empirical A–Q curve (rectangular hyperbola)
% ------------------------------------------------------------
para0 = [0.1 30 1 0.5];  
[para_fit_aq, ~] = lsqcurvefit(@A_Q_curve, para0, A_Q_am(:,1), A_Q_am(:,2), [0 0 0 0], [inf inf inf 1]);
PPFD_m = 0:1:2500;
for i=1:length(PPFD_m)
    A_m(i) = A_Q_curve(para_fit_aq, PPFD_m(i));
end
figure; plot(PPFD, A,'*'); hold on; plot(PPFD_m, A_m,'-');
xlabel('PPFD'); ylabel('A')

% ------------------------------------------------------------
% Record metadata (SPAD and genotype)
% ------------------------------------------------------------
spad_all = readtable(name, 'Sheet','spad');
spad = table2array(spad_all(1,1));
genotype_aq  = repmat(table2array(spad_all(1,2)), size(A_Q_am,1), 1);
genotype_aci = repmat(table2array(spad_all(1,2)), size(A_Ci_am,1), 1);
spad_aq  = repmat(spad, size(A_Q_am,1), 1);
spad_aci = repmat(spad, size(A_Ci_am,1), 1);

aq_table  = table(A_Q_am, spad_aq, genotype_aq);
aci_table = table(A_Ci_am, spad_aci, genotype_aci);

% Final parameter summary
result = [gbs, vcmax, vpmax, gm];
para_result = [result, R_light, para_fit_atp, gbs_jj/1000, para_fit_gs, para_fit_aq];
