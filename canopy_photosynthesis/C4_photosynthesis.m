function [A, para_temperature, A_all] = C4_photosynthesis(PPFD_in, Ci_in, Oi, Tleaf, para)
%% =========================================================================
% Function: C4_photosynthesis
% Purpose : Compute net CO2 assimilation rate (A) for a single leaf in C4 plants
% Inputs  :
%   PPFD_in - Photosynthetic photon flux density (µmol/m²/s)
%   Ci_in   - Intercellular CO2 concentration (ppm)
%   Oi      - Oxygen concentration (fraction)
%   Tleaf   - Leaf temperature (°C)
%   para    - Parameter vector: [Jatpmax25, YJatpll, theta_e, Rd25, Vcmax25, Vpmax25, gbs25, gm25]
% Outputs :
%   A                 - Net CO2 assimilation (µmol CO2/m²/s)
%   para_temperature  - Temperature-adjusted key parameters [Jatpmax, Vcmax, Vpmax, Rd, gbs, gm]
%   A_all             - All four potential assimilation rates [AEE, AET, ATT, ATE]
%% =========================================================================

%% Environment parameters
PPFD = PPFD_in.*10^(-6);
Ci =Ci_in.*10^(-6);
% Oi = para(3);
% Tleaf=para(4);
%% Electronic transport parameters
Jatpmax25 = para(1).*10^(-6); %  light-saturated ATP production rate (mole m-2 s-1)
YJatpll = para(2);% initial(or maximum) quantum yield for ATP production(i.e. conversion efficieny of PPFD into Jatp)
theta_e = para(3); %Curvature of the nonrectangular hyperbola describing the PPFD dependence of JATP dimensinless 
Rd25 =para(4).*10^(-6);  % Day respiration at base temperature (RTEMP) (mol m-2 s-1)
%% Enzyme Parameters 
Vcmax25 = para(5)*10^(-6); %Maximum rate of Rubisco carboxylation (mol m-2 s-1) at 25C
Vpmax25 =para(6)*10^(-6);% The maximum rate of PEP carboxylation (mol m-2 s-1) at 25C

%% cell conductance parameters
gbs25=para(7)*10^(-3);%1.47*10^(-3); %2.87*10^(-3) # Bundle sheath conductance (mol m-2 s-1) at 25c 2016 jxb
gm25 =para(8);%1.78; %   1.33  mesophyll conductance (mol m-2 s-1) at 25c 2016 jxb

%% 
uoc25 =0.047;  % the coefficient that lumps diffusivities and solubilities of CO2 and O2 in water 2016 JXB

%% electronic transport
alpha_e =0.15; %0.1;  % Fraction of PSII active in bundle sheath dimensionless 2016 JXB
x = 0.4; % fraction of ATP allocated to the C4 cycle 2016 JXB
Gamma_e25 =0.000193; %0.0001747; %half the reciprocal of Rubisco specificity at 25C 2016 JXB

%% michaelis-menten constant 
kp25 = 80*10^(-6);%60*10^(-6);%40 *10^(-6) # michaelis-menten constant of PEPc for co2 at 25c (bar) 2016 JXB 
kmc25 =650*10^(-6);%960*10^(-6);%485*10^(-6) # michaelis-menten constant of Rubisco for co2 at 25c (bar) 2016 JXB
kmo25 =450000*10^(-6);%497000*10^(-6);% #146000*10^(-6) # michaelis-menten constant of Rubisco for o2 at 25c (bar) 2016 JXB

%% activation and deactivation energy 2016 jxb
Evcmax = 53.4; %activation energy for Vcmax (kj mol-1)
Egamma = 27.4; % activation energy for half the reciprocal of Rubisco specificity (kj mol-1)
Ekmc =35.6;   %activation energy for kmc  (kj mol-1)
Ekmo = 15.1; %activation energy for kmo  (kj mol-1)
Evpmax = 37; %activation energy for vpmax  (kj mol-1)
Dvpmax = 214.5; %deactivation energy for vpmax  (kj mol-1)
Svpmax =0.663; %Entropy term for vpmax  (kj k-1 mol-1)
Ekp =68.1; %activation energy for kp (kj mol-1)
Egm =49.6; %activation energy for gm (kj mol-1)
Dgm = 437.4; %deactivation energy for gm (kj mol-1)
Sgm= 1.4; %Entropy term for gm (kj k-1 mol-1)
Egbs=116.77; %activation energy for gbs (kj mol-1)
Dgbs= 264.60; %deactivation energy for gbs (kj mol-1)
Sgbs = 0.86; %Entropy term for gbs (kj k-1 mol-1)
Erd= 41.85; %activation energy for Rd (kj mol-1)
Rugs = 0.008314; % universal gas constant (kj K-1 mol-1)

Ejmax = 32.8; %activation energy for vpmax  (kj mol-1)
Djmax = 220; %deactivation energy for vpmax  (kj mol-1)
Sjmax =0.703; %Entropy term for vpmax  (kj k-1 mol-1)
%%  Parameters modify by leaf temperature 
Vcmax = para_temperature_1(Vcmax25,Evcmax,Rugs,Tleaf) ;
Gamma_e= para_temperature_1(Gamma_e25,Egamma,Rugs,Tleaf); 
kp =para_temperature_1(kp25,Ekp,Rugs,Tleaf);
kmc =para_temperature_1(kmc25,Ekmc,Rugs,Tleaf);
kmo = para_temperature_1(kmo25,Ekmo,Rugs,Tleaf);
Rd = para_temperature_1(Rd25,Erd,Rugs,Tleaf);
Rm =0.5*Rd;

Vpmax = para_temperature_2(Vpmax25,Evpmax,Dvpmax,Svpmax,Rugs,Tleaf);
Jatpmax = para_temperature_2(Jatpmax25,Ejmax,Djmax,Sjmax,Rugs,Tleaf);
gm = para_temperature_2(gm25,Egm,Dgm,Sgm,Rugs,Tleaf);
gbs= para_temperature_2(gbs25,Egbs,Dgbs,Sgbs,Rugs,Tleaf);
para_temperature=[Jatpmax/10^(-6),Vcmax/10^(-6),Vpmax/10^(-6),Rd/10^(-6),gbs/10^(-3),gm25];
uoc= uoc25*exp(-(1.63/Rugs*(1/298-1/(273+Tleaf))));
%%  Functions to modify leaf nitrogen content
% gm <- gm_t*((-0.0007+1.989*leaf_N)/(-0.0007+1.989*1.1))
% gbs <- gbs_t*((1049.4-616.37*1.1)/(1049.4-616.37*leaf_N))

%% JATP
Jatp =(YJatpll.*PPFD+Jatpmax-sqrt((YJatpll.*PPFD+Jatpmax).^2-4.*theta_e.*Jatpmax.*YJatpll.*PPFD))./(2.*theta_e);

%% Parameters for enzyme(Rubisco)-limited rate and e- transport-limited rate
x1_enzyme = Vcmax;
x2_enzyme =kmc/kmo;
x3_enzyme = kmc;

x1_e = (1-x).*Jatp./3;
x2_e = 7.*Gamma_e./3;
x3_e = 0;
Vp_t = x.*Jatp./2;

%% ATE is the rate when the C4 cycle is limited by e- transport and the c3 cycle is limited by enzyme activity
ATE = quadratic_solution(Ci,Oi,x1_enzyme,x2_enzyme,x3_enzyme,Vp_t,uoc,alpha_e,Gamma_e,gm,gbs,Rm,Rd)*10^6;

%% ATT is the rate when both C4 and C3 cycles are limited by e- transport. 
ATT = quadratic_solution(Ci,Oi,x1_e,x2_e,x3_e,Vp_t,uoc,alpha_e,Gamma_e,gm,gbs,Rm,Rd)*10^6;

%% AEE is the net CO2 assimilation rate when both C4 and C3 cycles are limited by enzyme activity
AEE =cubic_solution(Ci,Oi,x1_enzyme,x2_enzyme,x3_enzyme,Vpmax,uoc,alpha_e,Gamma_e,gm,gbs,Rm,Rd,kp)*10^6;

%% AET is the net rate when the C4 cycle is limited by enzyme activity and the C3 cycle is limited by e- transport
AET = cubic_solution(Ci,Oi,x1_e,x2_e,x3_e,Vpmax,uoc,alpha_e,Gamma_e,gm,gbs,Rm,Rd,kp)*10^6;

%% the net CO2 assimilation rate
A = min([AEE,AET,ATT,ATE],[],2);
A_all=[AEE,AET,ATT,ATE];
end

%% Two functions to modify parameters by leaf temperature 
 function para_current=para_temperature_1(para_25,E,R,Tleaf) 
  para_current=para_25*exp(E/R*(1/298-1/(273+Tleaf)));
 end

function para_current =para_temperature_2(para_25,E,D,S,R,Tleaf)
  para_current=para_25*exp(E/R*(1/298-1/(273+Tleaf)))*(1+exp((S-D/298)/R))/(1+exp((S-D/(273+Tleaf))/R));
end

%% the quadratic solution to ATE or ATT
function A=quadratic_solution(Ci,Oi,x1,x2,x3,Vp,uoc,alpha_e,Gamma_e,gm,gbs,Rm,Rd)
a = x2.*gm.*alpha_e./uoc-gm-gbs;
b = (gm+gbs).*(x1-Rd)+gm.*(Ci.*gbs+Vp-Rm)+(x3+x2.*Oi).*gm.*gbs+(x1.*Gamma_e+x2.*Rd).*gm.*alpha_e/uoc;
c = -gm.*(Ci.*gbs+Vp-Rm).*(x1-Rd)+gm.*gbs.*(x1.*Gamma_e.*Oi+Rd.*(x3+x2.*Oi));
A=(-b +sqrt(b.^2-4.*a.*c))./(2.*a);
end

%% the cubic expression and its solution to AEE or AET
function A=cubic_solution(Ci,Oi,x1,x2,x3,Vpmax,uoc,alpha_e,Gamma_e,gm,gbs,Rm,Rd,kp)
  d = gm.*(Rm-Vpmax-Ci.*(gm+2.*gbs)-kp.*(gm+gbs));
  f = gm.^2.*(Ci.*Vpmax+(Ci+kp).*(gbs.*Ci-Rm));
  k = gm.^2.*gbs.*(Ci+kp);
  m = d-(x3+x2.*Oi).*gm.*gbs+(Rd-x1).*(gm+gbs)-(x1.*Gamma_e.*gm+x2.*Rd.*gm-x2.*k./gbs).*alpha_e/uoc;
  n = f+(x3+x2.*Oi).*k+d.*(Rd-x1)-gm.*gbs.*(x1.*Gamma_e.*Oi+Rd.*(x3+x2.*Oi))+(x1.*Gamma_e+x2.*Rd).*k.*alpha_e./(uoc.*gbs);
  o = Rd.*(f+(x3+x2.*Oi).*k)-x1.*(f-k.*Gamma_e.*Oi);
  p = m./(gm+gbs-x2.*gm.*alpha_e./uoc);
  q = n./(gm+gbs-x2.*gm.*alpha_e./uoc);
  r = o./(gm+gbs-x2.*gm.*alpha_e./uoc);
  Q = (p.^2-3.*q)./9;
  U = (2.*p.^3-9.*p.*q+27.*r)./54;
  Y = acos(U./sqrt(Q.^3));
  A = -2.*sqrt(Q).*cos(Y./3)-p./3;  
end
