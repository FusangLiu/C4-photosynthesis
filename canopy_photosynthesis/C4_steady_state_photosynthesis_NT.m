function [A_net, gs, gv, T_leaf, E, Ci, para_temp] = C4_steady_state_photosynthesis_NT(PPFD, CO2, Oi, T_air, RH, para, g0, g1, water_correct)
%% =========================================================================
% Function: C4_steady_state_photosynthesis_NT
% Purpose : Compute C4 photosynthesis without dynamic leaf temperature (T_leaf = T_air)
% Inputs  :
%   PPFD        - Photosynthetic photon flux density (µmol/m²/s)
%   CO2         - Ambient CO2 concentration (µmol/mol)
%   Oi          - Ambient O2 concentration (mol/mol)
%   T_air       - Air temperature (°C)
%   RH          - Relative humidity (%)
%   para        - Photosynthesis parameters [Jmax, kll, theta, Rd, Vcmax, Vpmax, gbs, gm]
%   g0, g1     - Stomatal conductance parameters
%   water_correct - Correction factor for water availability (unused here)
% Outputs :
%   A_net       - Net photosynthetic rate (µmol/m²/s)
%   gs          - Stomatal conductance (mol/m²/s)
%   gv          - Total leaf conductance to water vapor (mol/m²/s)
%   T_leaf      - Leaf temperature (°C, set equal to T_air)
%   E           - Transpiration rate (mmol/m²/s)
%   Ci          - Intercellular CO2 concentration (µmol/mol)
%   para_temp   - Updated photosynthesis parameter values after temperature correction
%% =========================================================================

iters = 1;               % Iteration counter
T_leaf = T_air;          % Leaf temperature fixed at air temperature
Ci = 400 * 0.7;          % Initial guess for intercellular CO2
Ci_old = 400;            % Previous iteration Ci

%% ----------------- Photosynthesis computation for PPFD > 20 -----------------
if PPFD > 20
    % Iterate until convergence of Ci (intercellular CO2) or max iterations reached
    while abs(Ci_old - Ci) > 5 && iters < 100
        iters = iters + 1;
        
        % Compute photosynthesis rate at fixed T_leaf
        [A_net, para_temp] = C4_photosynthesis(PPFD, Ci, Oi, T_air, para);
        
        % Compute stomatal conductances using Ball-Berry model
        [gs, gb, gh] = stomatal_conductance_BB_est([g0, g1], [A_net, CO2, RH, T_air]);
        gs = max(gs, 1e-9);  % avoid zero conductance
        gb = max(gb, 1e-9);
        
        % Total leaf conductance to water vapor
        gv = gs * gb / (gs + gb);
        
        % Compute transpiration rate (leaf temp = air temp)
        [~, E] = energy_balance(PPFD, T_air, RH / 100, gv, gh);
        
        % Update Ci
        Ci_old = Ci;
        Ci = CO2 - A_net / gv;
        
        % Check for non-physical Ci values
        if Ci <= 0
            warning('Ci iteration dropped to non-physical value, using previous Ci=%.4f', Ci_old);
            Ci = Ci_old;
            break;
        end
    end

%% ----------------- Low light condition (PPFD <= 20) -----------------
else
    % Direct computation without iteration
    [A_net, para_temp] = C4_photosynthesis(PPFD, Ci, Oi, T_air, para);
    [gs, gb, gh] = stomatal_conductance_BB_est([g0, g1], [A_net, CO2, RH, T_air]);
    gs = g0;                 % Minimum stomatal conductance
    gv = gs * gb / (gs + gb); % Total leaf conductance
    
    % Compute transpiration rate
    [T_leaf, E] = energy_balance(PPFD, T_air, RH / 100, gv, gh);
end

% Ensure net photosynthesis is real-valued
A_net = real(A_net);

end
