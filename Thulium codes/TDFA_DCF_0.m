%close all
%clearvars
%clc
% addpath('D:/PhD/M-codes/Functions')

function [Pp_out,Z,Spectrum_out,Ps_out,lambda_s,Excited_state_pop_out,Ps,N1,N2,bp,bs] = TDFA_DCF_0(L,N0,NAs,NAp,neff,rs,rp,alpha,lambda_p,lambda_s_start,lambda_s_end,lambda_s_step,tau21,tau32,K1,Pump_power, Signal_power,inserting_data,dz)
%% Physical parameters

h=6.626e-34;        %planck's constant in Js
c=3e8;              %speed of light in meter per sec
AN=6.02e23;         %Avogadro no

%% TDF parameters

fiber.L = L; %2;                    % Length of gain fiber
fiber.N0 = N0; %2.74e26;             % Total dopant concentration in /m^3
fiber.NAs = NAs; %0.15;               % Numerical Aperture of the core
fiber.NAp = NAp; %0.46;               % Numerical Aperture of the inner clad
fiber.neff = neff; %1.5;               % refractive index
fiber.rs = rs; %5e-6;                % radius of core in meters
fiber.rp = rp; %65e-6;               % radius of clad in meters

fiber.alpha = alpha; %0;

%% Pump and signal channels

fiber.Lp = lambda_p*1e-9; %790e-9;

l1 = lambda_s_start; %1568e-9;
l2 = lambda_s_end; %2155e-9;
dl = lambda_s_step; %1e-9;
fiber.Ls = (l1:dl:l2)';

g = length(fiber.Ls);

[fiber.wp,~] = mfd(fiber.Lp,fiber.rp,fiber.NAp);
[fiber.ws,fiber.gammas] = mfd(fiber.Ls,fiber.rs,fiber.NAs);
fiber.gammap = (fiber.rs/fiber.rp)^2;

%% Spectroscopic parameters, absorption and emission cross sections

fiber.tau21 = tau21; %334.7e-6;         % Lifetime of level 2 in sec
fiber.tau32 = tau32; %14.2e-6;          % Lifetime of level 3 in sec

fiber.K1 = K1; %1.25e-22;            % K3212 cross relaxation coefficient
fiber.K2 = 0.084*fiber.K1;      % K2321 cross relaxation coefficient

sec = xlsread('Emission_crosssection_TDF.xls');    %emission cross-section at signal wavelength in square meter
sac = xlsread('Absorption_crosssection_TDF.xlsx');      %absorption cross-section at signal wavelength in square meter

le = sec(:,1)*1e6;
se = sec(:,2)*1e25;

la = sac(:,1)*1e6;
sa = sac(:,2)*1e25;

fiber.sap = 1e-25*fiber.gammap.*interp1(la,sa,fiber.Lp*1e6);

fiber.sas = 1e-25*fiber.gammas.*interp1(la,sa,fiber.Ls*1e6);
fiber.ses = 1e-25*fiber.gammas.*interp1(le,se,fiber.Ls*1e6);

fiber.PoA=2*h*(c^2)*dl./(fiber.Ls.^3);            %Local noise power

clear sec sac le se la sa

%% Initial conditions, array initialization

dz = dz; %1e-3;
z = 0:dz:fiber.L;

disp(size(fiber.Ls));
options = odeset('Maxstep',dz);

% Pump = 1e-3*input('Pump power (mW) = ');
Pump = Pump_power; %1900;
Signal = Signal_power; %00e-6;
sig_lam = 1900e-9;

S_in = zeros(1,g);
S_in(abs(fiber.Ls-sig_lam)<(dl/2)) = Signal;

[Z,Sol] = ode45(@amp3,z,[Pump S_in]',options,fiber);

Z=Z';

Pp = Sol(:,1)';
Ps = Sol(:,2:end)';
Spectrum = Ps(:,end)';              % Output spectrum
Po = sum(Spectrum);                 % Total output power (signal)
Power = sum(Ps);                    % Total signal power across fiber

%semilogy(fiber.Ls,Ps(:,end)')

%% Initial conditions for Q-switching

[N1,N2,bp,bs] = inv3(Pp,Ps,fiber);


lambda_s = fiber.Ls;

Ps_out = Power*1e3;
Pp_out = Pp*1e3;
Spectrum_out = 10*log10(Spectrum*1e3);
Excited_state_pop_out = [N1;N2]/fiber.N0;

if inserting_data == "y"
    
    disp(Pp);
    dlmwrite('Pp.csv',Pp);
    dlmwrite('Ps.csv',Ps);
    dlmwrite('N1.csv',N1);
    dlmwrite('N2.csv',N2);
    dlmwrite('bp.csv',bp);
    dlmwrite('bs.csv',bs);
    dlmwrite('lams.csv',fiber.Ls);

    disp("Initial conditions written into the .csv files - Open");

else 

    figure()
    subplot(221),plot(Z,Pp*1e3),xlabel('Fiber position (m)'),ylabel('Pump power (mW)')%,Z,Power*1e3),title('Power vs Length'),legend('Pump','Signal'),grid on
    subplot(222),plot(Z,Power*1e3),xlabel('Fiber position (m)'),ylabel('Signal power (mW)')
    subplot(223),plot(Z,[N1;N2]/fiber.N0),axis([0 fiber.L 0 1]),xlabel('Fiber position (m)'),ylabel('N_1/N_0, N_2/N_0')%,title('Excited State Population'),grid on
    subplot(224),plot(fiber.Ls,10*log10(Spectrum*1e3)),xlabel('Wavelength (m)'),ylabel('Power (dBm)')%,lam_s,Noise)
end



