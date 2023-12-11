%close all
%clearvars  
%clc

function [Pp_out,Z,Spectrum_out,Ps_out,lambda_s,Ps,N1,bp,bs,lams] = YDF_Amplifier_0(L,dz,N0,NAs,neff,rs,alpha,lambda_p,lambda_s_start,lambda_s_end,lambda_s_step,tauf,wdm_l_in,iso_l_in,Pp_in,Ps_in,file_path,inserting_data)

%% Excited_state_pop_out
%%addpath('D:\Q-switchedLaser_Codes')
%% Physical parameters

h = 6.6256e-34;        %planck's constant in Js
c = 3e8;              %speed of light in meter per sec
AN = 6.02e23;         %Avogadro no

%% EDF/YDF parameters

%inserting_data = "n";

fiber.L = L; %3.5;                  % Length of gain fiber
fiber.N0 = N0; %2.0e25;              % Total dopant concentration in /m^3
fiber.NAs =  NAs; %0.29; 0.11-->YDF               % Numerical Aperture of the core
fiber.NAp = 0.46;               % Numerical Aperture of the inner clad
fiber.neff = neff; %1.5;               % refractive index
fiber.rs = rs; %4.4e-6/2; 3e-6--->YDF           % radius of core in meters
fiber.rp = 65e-6;               % radius of clad in meters

fiber.alpha = alpha; %0;


%% Pump and signal channels

fiber.Lp = lambda_p*1e-9; %976e-9;            % Pump wavelength

l1 = lambda_s_start; % 1501e-9;  1000e-9--->YDF                 % Signal start wavelength
l2 = lambda_s_end; %1600e-9;   1100e-9--->YDF                % Signal stop wavelength
dl = lambda_s_step; %.1e-9;                      % Signal wavelength separation
fiber.Ls = (l1:dl:l2)';        % Signal wavelength array
lams = fiber.Ls;
%fiber.Ls = fiber.Ls * 1e-9;
%disp(['l1 = ',num2str(l1),' mW'])
g = length(fiber.Ls);           % Number of signal wavelength channels

[fiber.wp,fiber.gammap] = mfd(fiber.Lp,fiber.rs,fiber.NAs); % Mode-field radius and overlap factor for pump
[fiber.ws,fiber.gammas] = mfd(fiber.Ls,fiber.rs,fiber.NAs); % Mode-field radius and overlap factor for signal
% fiber.gammap = (fiber.rs/fiber.rp)^2;     % Required in case of clad-pumping

%% Spectroscopic parameters, absorption and emission cross sections

fiber.tau_f = tauf; %10e-3;            % Lifetime of metastable state

Data = dlmread(file_path); %'HG980.csv');    % Spectroscopic data of the gain fiber
%Data = Data';
fiber.sap = 1e-24*fiber.gammap.*Data(1,2);  % Pump absorption cross section   1e-25 for edfa
fiber.sep = 1e-24*fiber.gammap.*Data(1,3);  % Pump emission cross section


sa = Data(2:end,2);
se = Data(2:end,3);
ll = Data(2:end,1);

%data1 = Data(:,1);
%ll = data1(2:end);
fiber.sas = 1e-25*fiber.gammas.*interp1(ll,sa,fiber.Ls*1e9);    % Signal absorption cross section
fiber.ses = 1e-25*fiber.gammas.*interp1(ll,se,fiber.Ls*1e9);    % Signal emission cross section

fiber.PoA=2*h*(c^2)*dl./(fiber.Ls.^3);      % ASE noise power

clear sec sac ll se sa

%% components

wdm_p = wdm_l_in/4.343; %1.55/4.343;         % insertion loss of WDM at pump
iso_l = iso_l_in/4.343; %0.9/4.343;          % insertion loss of isolator
wdm_l = 1.1/4.343;          % insertion loss of WDM

%% Initial conditions, array initialization

%dz = 0.015; %7.5e-3;                    % Space step size
z = 0:dz:fiber.L;               % Space array over gain fiber length

options = odeset('Maxstep',dz);     % Step size setting for ODE solver

Pump = Pp_in*1e-3; %1e-3*input('Pump power (mW) = ');    % Pump power input
P_in = Pump*exp(-wdm_p);                    % Pump power entering the gain fiber
Signal = Ps_in*1e-3; %input('Signal power (mW) = ');% Signal power input (default value 0)
if(isempty(Signal))
    Signal = 0;
end
sig_lam = 1550e-9;                          % Signal wavelength

S_in = zeros(1,g);
u = find(abs(fiber.Ls-sig_lam)<(dl/2));
S_in(u) = Signal*exp(-wdm_l); % Initial conditions for all signal wavelengths

%% Solution

[Z,Sol] = ode45(@amp_2level,z,[P_in S_in],options,fiber);	% Solution of coupled rate and power propagation equations

Z=Z';       

Pp = Sol(:,1)';         % Solution of pump power over fiber length
Ps = Sol(:,2:end)';     % Solution of all signal wavelengths over fiber length

Spectrum = Ps(:,end)'*exp(-iso_l);  % Output spectrum
Po = sum(Spectrum);                 % Total output power (signal)
Power = sum(Ps);                    % Total signal power across fiber

%% Inversion and gain coefficients

[N1,bp,bs] = inv_2level(Pp,Ps,fiber);   % Ground-state population and gain coefficients of pump and signal
N2 = fiber.N0 - N1;                     % Metastable state population

gain = trapz(z,bs')*4.343;              % Gain spectrum (integrated over fiber length)
gs = 4.343*bs(u,:);                     % Gain coefficient at signal wavelength across the fiber

%Y = Pp;
%X = Z;
lambda_s = fiber.Ls;

Ps_out = Power*1e3;
Pp_out = Pp*1e3;
Spectrum_out = 10*log10(Spectrum*1e3);
Excited_state_pop_out = [N1;N2]/fiber.N0;

%% Output display and plots

%disp(['Pump coupled in EDF = ',num2str(P_in*1e3),' mW'])
%disp(['Signal wavelength = ',num2str(sig_lam*1e9),' nm'])
%disp(['Input Power = ',num2str(Signal*1e3),' mW'])
%disp(['Output Power = ',num2str(Po*1e3),' mW'])
%disp(['Residual Pump = ',num2str(Pp(end)*1e3),' mW'])
%disp(['Signal Gain = ',num2str(gain(u)),' dB'])
%disp(['Pump absorption = ',num2str((-10*log10((Pp(end))/P_in))/fiber.L),' dB/m'])

if inserting_data == "y"
    

    dlmwrite("Pp.csv",Pp*1e3);
    dlmwrite("Ps.csv",Ps*1e3);
    dlmwrite("N1.csv",N1);
    dlmwrite("bp.csv",bp);
    dlmwrite("bs.csv",bs);
    dlmwrite("lams.csv",fiber.Ls);

    disp("Initial conditions written into the .csv files - Open");
    
else
    figure()
    subplot(221),plot(Z,Pp*1e3),xlabel('Fiber position (m)'),ylabel('Pump power (mW)')%,Z,Power*1e3),title('Power vs Length'),legend('Pump','Signal'),grid on
    subplot(222),plot(Z,Power*1e3),xlabel('Fiber position (m)'),ylabel('Signal power (mW)')
    subplot(223),plot(Z,[N1;N2]/fiber.N0),axis([0 fiber.L 0 1]),xlabel('Fiber position (m)'),ylabel('N_1/N_0, N_2/N_0')%,title('Excited State Population'),grid on
    subplot(224),plot(fiber.Ls,10*log10(Spectrum*1e3)),xlabel('Wavelength (m)'),ylabel('Power (dBm)')%,lam_s,Noise)
end


%figure(),plot(fiber.Ls*1e9,gain),xlabel('Wavelength (m)'),ylabel('Gain (dB)')
%figure(),plot(Z,gs),xlabel('Fiber position (m)'),ylabel('Signal gain (dB)')

