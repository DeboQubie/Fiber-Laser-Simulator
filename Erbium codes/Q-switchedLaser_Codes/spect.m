function out = spect(lam,EDFA_N0,EDFA_w,NA)
%% initialization of physical constants
nf=1.5; % refractive index of EDFA
c=3e8;  % speed of light
v=c/nf; % velocity of light in EDFA
h=6.625e-34; % Planck's Constant
%% EDFA parameters

Data = dlmread('HG980.csv');
% lam = dlmread('lam.csv');

EDFA_lam_p = Data(1,1)*1e-9;
sig_a = Data(1,2)*1e-25;
sig_e = Data(1,3)*1e-25;

% EDFA_lam_p = lam(1)*1e-9;                          % Pump wavelength
% pump_abs = 250/4.343;                           % absorption coefficient for pump
% sig_a = interp1(Data(:,1),Data(:,2),EDFA_lam_p*1e9)*1e-25;          % absorption cross section for pump
% sig_e = interp1(Data(:,1),Data(:,3),EDFA_lam_p*1e9)*1e-25;          % emission cross section for pump
% EDFA_lam_p = 0.915e-6;                          % Pump wavelength
% pump_abs = 85/4.343;
% sig_a = 0.57e-24;
% sig_e = 0.019e-24;

[EDFA_wp,gamma_p] = mfd(EDFA_lam_p,EDFA_w,NA);
% gamma_p = 0.5;
% EDFA_N0 = 3.61e25;%pump_abs/sig_a;                               % Doping concentration
EDFA_ala_p = sig_a*EDFA_N0*gamma_p;                   % abs.coeff at pump
EDFA_ale_p = sig_e*EDFA_N0*gamma_p;

lam_s = lam(2:end);
dl = lam_s(2) - lam_s(1);
sig_as = interp1(Data(2:end,1),Data(2:end,2),lam_s);
sig_es = interp1(Data(2:end,1),Data(2:end,3),lam_s);

EDFA_lam_s = lam_s'*1e-9;                            % signal wavelengths
[EDFA_ws,gamma_s] = mfd(EDFA_lam_s,EDFA_w,NA);
% gamma_s(:,:) = 0.5;
EDFA_ala_s = sig_as'*1e-25*EDFA_N0.*gamma_s;                % abs. coeffs at signal wavelengths
EDFA_ale_s = sig_es'*1e-25*EDFA_N0.*gamma_s;                % emission coeffs at signal wavelengths

EDFA_P0 = h*c^2./EDFA_lam_s.^3*(dl*1e-9) * 2;           % spectral emission

t_f = 10e-3;                                     % flourescence lifetime of upper lasing level
out(:,1) = 1e-9*lam';
out(:,2) = [EDFA_wp EDFA_ws']';
out(:,3) = [EDFA_ala_p EDFA_ala_s']';
out(:,4) = [EDFA_ale_p EDFA_ale_s']';
out(:,5) = [t_f EDFA_P0']';
% p = length(EDFA_lam_s);
clear Data;