%close all
%clearvars
%clc

function [Z,lambda_s,Spectrum_out,Gain_out,Pp_out,Ps_out,fiber] = YDF_CWLaser_0(L,N0,NAs,neff,rs,alpha,lambda_p,lambda_s_start,lambda_s_stop,lambda_s_step,tauf,wdm_p_in,iso_l_in,cou_l_in,c_out_in,wdm_l_in,aom_l_in,c_l_in,Pp,file_name,inserting_data)

%addpath('D:/PhD/M-codes/Functions')

%% Physical parameters

h=6.626e-34;        %planck's constant in Js
c=3e8;              %speed of light in meter per sec
AN=6.02e23;         %Avogadro no

%% EDF/YDF parameters

fiber.L = L; %3.5;                  % Length of gain fiber
fiber.N0 = N0; %2.0e25;              % Total dopant concentration in /m^3
fiber.NAs = NAs; %0.29;               % Numerical Aperture of the core
fiber.NAp = 0.46;               % Numerical Aperture of the inner clad
fiber.neff = neff; %1.5;               % refractive index
fiber.rs = rs; %4.4e-6/2;            % radius of core in meters
fiber.rp = 65e-6;               % radius of clad in meters

fiber.alpha = alpha; %0;


%% Pump and signal channels

fiber.Lp = lambda_p*1e-9; %976e-9;            % Pump wavelength

l1 = lambda_s_start; %1501e-9;  %1501e-9;                   % Signal start wavelength
l2 = lambda_s_stop; %1600e-9;  %1600e-9;                   % Signal stop wavelength
dl = lambda_s_step; %0.5e-9;                      % Signal wavelength separation
fiber.Ls = (l1:dl:l2)';         % Signal wavelength array


g = length(fiber.Ls);           % Number of signal wavelength channels

[fiber.wp,fiber.gammap] = mfd(fiber.Lp,fiber.rs,fiber.NAs); % Mode-field radius and overlap factor for pump
[fiber.ws,fiber.gammas] = mfd(fiber.Ls,fiber.rs,fiber.NAs); % Mode-field radius and overlap factor for signal
% fiber.gammap = (fiber.rs/fiber.rp)^2;     % Required in case of clad-pumping

%% Spectroscopic parameters, absorption and emission cross sections

fiber.tau_f = tauf; %10e-3;            % Lifetime of metastable state

Data = dlmread(file_name); %'HG980.csv');    % Spectroscopic data of the gain fiber

fiber.sap = 1e-25*fiber.gammap.*Data(1,2);  % Pump absorption cross section
fiber.sep = 1e-25*fiber.gammap.*Data(1,3);  % Pump emission cross section

ll = Data(2:end,1);
sa = Data(2:end,2);
se = Data(2:end,3);

fiber.sas = 1e-25*fiber.gammas.*interp1(ll,sa,fiber.Ls*1e9);    % Signal absorption cross section
fiber.ses = 1e-25*fiber.gammas.*interp1(ll,se,fiber.Ls*1e9);    % Signal emission cross section

fiber.PoA=2*h*(c^2)*dl./(fiber.Ls.^3);      % ASE noise power

clear sec sac ll se sa

%% components

wdm_p = wdm_p_in; %1.55;         % Insertion loss of WDM at pump
iso_l = iso_l_in; %0.9;          % Insertion loss of isolator
cou_l = cou_l_in; %3.35;         % Insertion loss of coupler at feedback port
c_out = c_out_in; %3.35;         % Insertion loss of coupler at output port
wdm_l = wdm_l_in; %1.1;          % Insertion loss of WDM
aom_l = aom_l_in; %0;            % Insertion loss of intracavity modulator : NOT REQUIRED          
c_l = c_l_in; %0;              % Additional connector loss (use it to vary the cavity loss)

Cav_Loss = iso_l + cou_l + wdm_l + aom_l + c_l;     % Total cavity Loss in dB
loss = exp(-Cav_Loss/4.343);                        % Cavity loss factor 
rout = exp(-(iso_l+c_out+c_l)/4.343);               % Output coupling factor


%% Initial conditions, array initialization

dz = 0.75e-2;                  % Space step size
z = 0:dz:fiber.L;           % Space array over gain fiber length

options = odeset('Maxstep',dz);     % Step size setting for ODE solver

Pump = Pp*1e-3; %1e-3*input('Pump power (mW) = ');    % Pump power input
P_in = Pump*exp(-wdm_p/4.343);                    % Pump power entering the gain fiber
S_in = zeros(1,g);                          % Signal power initalized as zero


%% Laser loop

eps = 1e-2;            % Relative error tolerance = 1%
err = 1;                % Initial error value (initalized above tolerance)
r=0;                    % Loop counter
Pp = 0;Ps = 0;          % Arbitrary initial values of pump and signal powers
Po = zeros(1,10000);    % Array to store output power values
tic
while(err>eps)
    r = r + 1;          % Counter incremented
    PP = Pp; PS = Ps;   % Pump and signal powers from previous round-trip copied, to compare with the powers evaluated in current round-trip
    
    [Z,Sol] = ode45(@amp_2level,z,[P_in S_in],options,fiber);	% Solution of coupled rate and power propagation equations
    
    Z=Z';
    
    Pp = Sol(:,1)';         % Solution of pump power over fiber length
    Ps = Sol(:,2:end)';     % Solution of all signal wavelengths over fiber length
    
    Power = sum(Ps);        % Total signal power across fiber
    P_end = Ps(:,end)';     % Signal power (spectrum) at the end of gain fiber
    Spectrum = rout*P_end;  % Spectrum of the output of laser cavity
    S_in = loss*P_end;      % Spectrum of signal after undergoing cavity loss : to be provided as initial value for next round-trip
    
    [N1,bp,bs] = inv_2level(Pp,Ps,fiber);    % Ground-state population and gain coefficients of pump and signal
    
    err = max([max(abs(Pp-PP)./PP),max(max(abs(Ps-PS)./PS))]);  % Relative error between pump and signal powers in current and previous round-trips
    
    disp(r)
%     fprintf('Round-trip : %u, relative error = %2.2f %%\n',r,err*100)   % Displaying round-trip count and relative error
    
    Po(r) = sum(Spectrum);          % Total output power
    
    gain = trapz(z,bs')*4.343;      % Gain spectrum (integrated over fiber length)
    
    subplot(221),plot(Z,Pp*1e3)%,xlabel('Length (m)'),ylabel('Power (mW)')%,Z,Power*1e3),title('Power vs Length'),legend('Pump','Signal'),grid on
    subplot(222),plot(Z,Power*1e3)%,xlabel('Length (m)'),ylabel('Power (mW)')
    subplot(223),plot(fiber.Ls*1e9,gain,fiber.Ls*1e9,Cav_Loss*ones(1,g))%,xlabel('Wavelength (nm)'),ylabel('Gain (dB)')
    subplot(224),plot(fiber.Ls*1e9,10*log10(Spectrum*1e3))%,xlabel('Wavelength (m)'),ylabel('Power (dBm)')%,lam_s,Noise)
    pause(5e-3)

    lambda_s = fiber.Ls*1e9;
    Spectrum_out = 10*log10(Spectrum*1e3);
    Gain_out = gain;
    Pp_out = Pp*1e3;
    Ps_out = Power*1e3;

end
toc
%% Displaying results

if inserting_data == "y"
    dlmwrite('Pp.csv',Pp*1e3);
    dlmwrite('Ps.csv',Ps*1e3);
    dlmwrite('N1.csv',N1);
    dlmwrite('bp.csv',bp);
    dlmwrite('bs.csv',bs);
    dlmwrite('lams.csv',fiber.Ls);

    disp("Initial conditions written into the .csv files - Closed");

else
    disp(['CW Power = ',num2str(Po(r)*1e3),' mW = ',num2str(10*log10(Po(r)*1e3)),' dBm'])
    disp(['Lasing wavelength = ',num2str(fiber.Ls(Spectrum==max(Spectrum))*1e9),' nm'])

end


    