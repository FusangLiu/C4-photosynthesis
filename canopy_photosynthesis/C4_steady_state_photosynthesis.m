function [A_net, gs, gv, T_leaf, E, Ci, para_temp] = ...
    C4_steady_state_photosynthesis(PPFD, CO2, Oi, T_air, RH, para, g0, g1, water_correct)
%% =========================================================================
% Function: C4_steady_state_photosynthesis
% Purpose : Calculate steady-state C4 photosynthesis for a single leaf facet
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
% Inputs  :
%   PPFD          - Photosynthetic photon flux density (µmol/m²/s)
%   CO2           - Ambient CO2 concentration (ppm)
%   Oi            - O2 concentration (fraction, e.g., 0.21)
%   T_air         - Air temperature (°C)
%   RH            - Relative humidity (%)
%   para          - Photosynthesis parameters [Jmax, kLL, theta_e, Rd, Vcmax, Vpmax, gbs, gm]
%   g0, g1        - Ball-Berry stomatal conductance parameters
%   water_correct - Water stress factor (1 = no stress)
% Outputs :
%   A_net     - Net photosynthesis (µmol CO2/m²/s)
%   gs        - Stomatal conductance (mol/m²/s)
%   gv        - Total leaf conductance to water vapor (mol/m²/s)
%   T_leaf    - Leaf temperature (°C)
%   E         - Transpiration rate (mol/m²/s)
%   Ci        - Intercellular CO2 concentration (ppm)
%   para_temp - Photosynthesis intermediate parameters
%% =========================================================================

iters = 1;            % Iteration counter
T_leaf = T_air + 1;   % Initial guess for leaf temperature
T_leaf_old = T_air;   % Store previous leaf temperature for convergence check

Ci = 400 * 0.7;       % Initial guess of intercellular CO2 concentration

%% ---------------- Step 1: High light scenario ----------------
if PPFD > 30
    % Iterate until leaf temperature and photosynthesis converge
    while abs(T_leaf_old - T_leaf) > 0.1 && iters < 100
        iters = iters + 1;
        T_leaf_old = T_leaf;
        
        % Compute leaf photosynthesis using C4 model
        [A_net, para_temp] = C4_photosynthesis(PPFD, Ci, Oi, T_leaf, para);
        
        % Estimate stomatal conductance using Ball-Berry type model
        [gs, gb, gh] = stomatal_conductance_BB_est([g0, g1], [A_net, CO2, RH, T_leaf]);
        
        % Avoid nonphysical negative conductances
        gs = max(gs, 1e-9);
        gb = max(gb, 1e-9);
        
        % Total leaf conductance to water vapor
        gv = gs * gb / (gs + gb);
        
        % Update leaf temperature using energy balance
        [T_leaf, E] = energy_balance(PPFD, T_air, RH/100, gv, gh);
        
        % Update intercellular CO2 concentration based on assimilation
        Ci_old = Ci;  % store previous value
        Ci = CO2 - A_net / gv;
        
        if Ci <= 0
            warning('Ci dropped to non-physical value, using previous Ci = %.4f', Ci_old);
            Ci = Ci_old;
            break;
        end
    end

%% ---------------- Step 2: Low light scenario ----------------
else
    % Low light (PPFD ≤ 30), photosynthesis limited
    [A_net, para_temp] = C4_photosynthesis(PPFD, Ci, Oi, T_air, para);
    
    % Stomatal conductance at low light
    [gs, gb, gh] = stomatal_conductance_BB_est([g0, g1], [A_net, CO2, RH, T_air]);
    
    % Set stomatal conductance to g0 minimum
    gs = g0;
    gv = gs * gb / (gs + gb);
    
    % Update leaf temperature and transpiration
    [T_leaf, E] = energy_balance(PPFD, T_air, RH/100, gv, gh);
end

%% ---------------- Step 3: Ensure real output ----------------
A_net = real(A_net);  % Remove any imaginary parts due to numeric instability

end
