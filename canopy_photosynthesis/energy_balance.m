function [t_l, E] = energy_balance(PPFD, t, RH, gv, gh)
%% =========================================================================
% Function: energy_balance
% Purpose : Compute leaf temperature and transpiration rate based on energy balance
% Inputs  :
%   PPFD - Photosynthetic photon flux density (µmol/m²/s)
%   t    - Air temperature (°C)
%   RH   - Relative humidity (%)
%   gv   - Total leaf conductance to water vapor (mol/m²/s)
%   gh   - Boundary layer conductance for heat (m/s or mol/m²/s)
% Outputs :
%   t_l  - Leaf temperature (°C)
%   E    - Transpiration rate (mmol H2O/m²/s)
%% =========================================================================

% ---------------- Step 1: Vapor pressure calculation ----------------
es = 0.611 * exp((17.502 * t) / (240.97 + t));  % Saturated vapor pressure at air temp, kPa
ea = es * RH;                                   % Ambient vapor pressure of air, kPa

% Adjust thermal conductance coefficient
gh = gh * 0.4;   % empirical adjustment for boundary layer heat transfer

% ---------------- Step 2: Solve leaf temperature ----------------
% Define function for energy balance residual (should be zero at equilibrium)
fun_t = @(t_l) est_leaf_temperature(t_l, t, gv, gh, ea, PPFD);

% Solve leaf temperature using bisection method over +/-30°C range
t_l = bisection(fun_t, t - 30, t + 30);

% ---------------- Step 3: Compute transpiration ----------------
es_l = 0.611 * exp((17.502 * t_l) / (240.97 + t_l)); % Saturated vapor pressure at leaf temp
delta_w = es_l - ea;                                 % Vapor pressure deficit, kPa
E = gv * delta_w;                                    % Transpiration, mmol/m²/s

end

%% =========================================================================
function res = est_leaf_temperature(t_l, t, gv, gh, ea, PPFD)
%% Purpose: Compute residual of energy balance at given leaf temperature
% Residual = Net radiation absorbed - Sensible heat - Latent heat
%% ---------------- Step 1: Constants ----------------
eps = 0.95;        % Leaf thermal emissivity
alpha_s = 0.85;    % Leaf absorption coefficient for PAR
sigma = 5.670e-8;  % Stefan-Boltzmann constant (W/m²/K⁴)
Lambda = 44;       % Latent heat of vaporization (kJ/mol)
Cp = 29.3;         % Specific heat of air (J/mol/K)
k = 1 / 4.55;      % Conversion factor from µmol/m²/s to W/m²

%% ---------------- Step 2: Temperature conversions ----------------
Tk_air = t + 240.97;   % Absolute air temperature in K
Tk = t_l + 240.97;     % Absolute leaf temperature in K
delta_T = t_l - t;     % Temperature difference leaf-air

%% ---------------- Step 3: Vapor pressure ----------------
es_l = 0.611 * exp((17.502 * t_l) / (240.97 + t_l)); % Saturated vapor pressure at leaf temp
PAR = PPFD * k;                                       % Convert PPFD to W/m²

%% ---------------- Step 4: Radiation components ----------------
R_sw = PAR * alpha_s;               % Absorbed shortwave radiation
R_wall = 2 * eps * sigma * Tk_air^4; % Thermal radiation absorbed from surroundings
R_leaf = 2 * eps * sigma * Tk^4;     % Thermal radiation emitted by leaf
R_thermal = R_wall - R_leaf;         % Net thermal radiation absorbed
R_net = R_sw + R_thermal;            % Total net radiation absorbed

%% ---------------- Step 5: Heat fluxes ----------------
delta_w = es_l - ea;       % Vapor pressure gradient
E = gv * delta_w;          % Latent heat flux (transpiration)
H = Cp * gh * delta_T;     % Sensible heat flux
Lambda_E = Lambda * E;     % Latent heat flux in W/m²

%% ---------------- Step 6: Energy balance residual ----------------
res = R_net - H - Lambda_E; % Residual: zero at equilibrium leaf temperature
end
