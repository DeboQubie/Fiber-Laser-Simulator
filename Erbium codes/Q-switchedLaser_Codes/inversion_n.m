function [N1,bp,bs] = inversion_n(Pp,Ps,Input,fiber)
%% initialization of physical constants

% nf=1.5; % refractive index of YDFA
c=3e8;  % speed of light
% v=c/nf; % velocity of light in YDFA
h=6.625e-34; % Planck's Constant
%% YDFA parameters

% Input = dlmread('input.csv');
% fiber = dlmread('fiber.csv');
YDFA_w = fiber(2);                                  % mode field radius
YDFA_lam_p = Input(1,1);                          % Pump wavelength

YDFA_wp = Input(1,2);

YDFA_N0 = fiber(1);%pump_abs/sig_a;                               % Doping concentration
YDFA_ala_p = Input(1,3);%sig_a*YDFA_N0*gamma_p;                   % abs.coeff at pump
YDFA_ale_p = Input(1,4);%sig_e*YDFA_N0*gamma_p;


YDFA_lam_s = Input(2:end,1);%lam_s'*1e-9;                            % signal wavelengths

YDFA_ws = Input(2:end,2);
YDFA_ala_s = Input(2:end,3);%sig_as'*1e-24*YDFA_N0.*gamma_s;                % abs. coeffs at signal wavelengths
YDFA_ale_s = Input(2:end,4);%sig_es'*1e-24*YDFA_N0.*gamma_s;                % emission coeffs at signal wavelengths

YDFA_P0 = Input(2:end,5);%h*c^2./YDFA_lam_s.^3*(dl*1e-9) * 2;           % spectral emission

t_f = Input(1,5);%;                                     % flourescence lifetime of upper lasing level

p = length(YDFA_lam_s);
clear Input fiber;
%% array initialization

r13 = (YDFA_ala_p/YDFA_N0)/(h*c/YDFA_lam_p)/(pi*YDFA_w^2);      % pump absorption rate
r31 = (YDFA_ale_p/YDFA_N0)/(h*c/YDFA_lam_p)/(pi*YDFA_w^2);      % pump emission rate
r12 = (YDFA_ala_s./YDFA_N0)./(h*c./YDFA_lam_s)./(pi*YDFA_w.^2);      % signal absorption rate
r21 = (YDFA_ale_s./YDFA_N0)./(h*c./YDFA_lam_s)./(pi*YDFA_w.^2);      % stimulated emission rate
A21 = 1/t_f;                                                    % spontaneous emission rate

%% function declaration

N1 = (r31*Pp + r21'*Ps + A21)./(r13*Pp + r31*Pp + r12'*Ps + r21'*Ps + A21);

bp = YDFA_ale_p*(1-N1) - YDFA_ala_p*N1;
bs = YDFA_ale_s*(1-N1) - YDFA_ala_s*N1;
