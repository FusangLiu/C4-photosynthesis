function [gs, gb, gh] = stomatal_conductance_BB_est(para, input)
%% =========================================================================
% Function: stomatal_conductance_BB_est
% Purpose : Estimate stomatal conductance using Ball-Berry model
% Inputs  :
%   para  - [g0, g1] stomatal parameters
%   input - [A_net, CO2, RH, T_air] for a leaf
% Outputs :
%   gs    - stomatal conductance to water vapor (mol/m²/s)
%   gb    - boundary layer conductance to water vapor (mol/m²/s/bar)
%   gh    - boundary layer conductance for heat (mmol/m²/s)
%% =========================================================================

g0 = para(1);  % Minimum stomatal conductance (mol/m²/s)
g1 = para(2);  % Sensitivity coefficient for Ball-Berry model

A_net = input(:,1); % Net photosynthesis rate (µmol/m²/s)
CO2   = input(:,2); % Ambient CO2 concentration (µmol/mol)
RH    = input(:,3) ./ 100; % Relative humidity fraction
T_air = input(:,4); % Air temperature (°C)

%% Compute boundary layer conductances
[gb, gbc, gh] = boundary_layer_conductance(T_air);

%% CO2 concentration at leaf surface (accounting for boundary layer resistance)
Cs = CO2 - A_net ./ gbc;

%% Ball-Berry stomatal model (simplified using RH at leaf surface)
% gs = g0 + g1 * (A_net * hs / Cs)
hs = RH; % assume leaf surface relative humidity equals ambient
gs = g0 + g1 .* (A_net .* hs ./ Cs);

end

%% =========================================================================
% Helper function to estimate relative humidity at leaf surface (unused here)
function res = est_hs(hs, g0, g1, gb, A_net, Cs, RH)
    gs = g0 + g1 .* (A_net .* hs ./ Cs);
    res = (hs - RH) .* gb - (1 - hs) .* gs;
end

%% =========================================================================
% Compute boundary layer conductances for water, CO2, and heat
function [gb, gbc, gh] = boundary_layer_conductance(T)
% Inputs: leaf temperature T (°C)
% Outputs:
%   gb  - boundary layer conductance to water (mol/m²/s/bar)
%   gbc - boundary layer conductance to CO2 (mol/m²/s/bar)
%   gh  - boundary layer conductance for heat (mmol/m²/s)

Tk_air = 240.97 + T; % Absolute temperature in Kelvin
u  = 2;              % Wind speed (m/s)
p  = 100;            % Air pressure (kPa)
w  = 10;             % Leaf width (cm)
sr = 1;              % Stomatal ratio
scr = (sr + 1)^2 / (sr^2 + 1); % Side conductance ratio
ocr = 1.4;           % Outdoor conductance ratio

%% Diffusion coefficients (m²/s)
Dw = 24.2e-6; % water vapor in air
Dc = 14.7e-6; % CO2 in air
Dh = 21.5e-6; % heat in air
Dm = 15.1e-6; % momentum in air

%% Convective heat transfer (Nusselt number)
d = 0.72 * w / 100; % characteristic leaf dimension (m)
Re = u * d / Dm;    % Reynolds number
Nu = 0.60 * sqrt(Re); % Nusselt number
gH = (Dh * Nu / d) * scr * ocr; % Convective heat conductance (m/s)

%% Convert to molar units
R = 8.314;          % J K^-1 mol^-1
gh = gH .* p ./ (R .* Tk_air) * 1000; % mmol/m²/s

%% Boundary layer conductance for water and CO2
rhw = (Dw / Dh)^(2/3);
gb = rhw .* gh / (p / 100);          % mol/m²/s/bar

drb = (Dw / Dc)^(2/3);
rbc = drb ./ gb;                     % resistance for CO2
gbc = 1 ./ rbc;                      % conductance for CO2 (mol/m²/s/bar)
end
