%clear all
%close all
%clc

function [Z_out,z_out,Pp_out,Ps_out,Pout_out,t_out,lambda_s_out,inv_out] = FTBS_CW_0(EDFA_N0,EDFA_L,EDFA_w,NA,lambda_p,tauf,L,iso_l_in,cou_l_in,aom_l_in,wdm_l_in,Ppf,Ppb,T,q,R,N,file_path)

%% initialization of physical constants

nf=1.5; % refractive index of EDFA
c=3e8;  % speed of light
v=c/nf; % velocity of light in EDFA
h=6.625e-34; % Planck's Constant
%% EDFA parameters
EDFA_N0 = EDFA_N0; %2.0e25;                               % Doping concentration
EDFA_L = EDFA_L; %2;                                    % length of EDFA
EDFA_w = EDFA_w; %4.4e-6/2;                               % EDFA Mode field radius
NA = NA; % 0.29;

Data = dlmread(file_path); %'HG980.csv');

EDFA_lam_p = lambda_p*1e-9; %0.98e-6;                            % Pump wavelength
[EDFA_wp,gamma_p] = mfd(EDFA_lam_p,EDFA_w,NA);
EDFA_ala_p = 1.55e-25*EDFA_N0*gamma_p;                   % abs.coeff at 980

EDFA_lam_s = Data(:,1)*1e-9;                  % signal wavelengths 1530nm,1540nm
[EDFA_ws,gamma_s] = mfd(EDFA_lam_s,EDFA_w,NA);

EDFA_ala_s = gamma_s.*Data(:,2)*EDFA_N0*1e-25;  % abs. coeffs at 1530,1540
EDFA_ale_s = gamma_s.*Data(:,3)*EDFA_N0*1e-25;  % emission coeffs at 1530,1540
EDFA_d_lam = EDFA_lam_s(2)-EDFA_lam_s(1);        % Noise photon bandwidth
EDFA_P0 = h*c^2./EDFA_lam_s.^3*EDFA_d_lam*2;   % spectral emission

t_f = tauf; %10e-3;                                     % flourescence lifetime of upper lasing level
p = length(EDFA_lam_s);

EDFA_A = pi*EDFA_w^2;                            % EDFA cross-section area
%% Cavity parameters

L = L; %10.0;                                        % length of the cavity
l1 = (L-EDFA_L)/5;                              % length of the fibre section between 2 components
tr = L/v;                                       % cavity roundtrip time
iso_l = iso_l_in/4.343; %0.9/4.343;                             % insertion loss of isolator
cou_l = cou_l_in/4.343; %3.3/4.343;                                % insertion loss of 3-dB coupler
aom_l = aom_l_in/4.343; %0.0/4.343;                                % insertion loss of AOM
% tsw = 150e-9;                                    % AOM switching time
wdm_l = wdm_l_in/4.343; %1.6/4.343;                             % insertion loss of WDM

Loss = (iso_l + cou_l + 25*aom_l + wdm_l)/L;
%% System Parameters

Pumpf = Ppf*1e-3; %0.600;                                   % Pump power
Pumpb = Ppb*1e-3; %0.0;

T = T; %1e-4;

% Rr = 10e3;
% T = 1/Rr;
% Th = 1e-6;
% Tl = T-Th;
q = q; %15;
%% Numerical parameters

R = R; %0.9;        % Stability factor : Courant parameter
N = N; %50;         % No. of space steps over EDFA
dz = EDFA_L/N;  % size of space steps
dt = R * dz/v;  % size of time steps
Nt = ceil(T/dt);      % No. of time steps
% Ntl = ceil(Tl/dt);
% Nth = ceil(Th/dt);

Nz = ceil(L/dz);      % No. of space steps over the cavity
nz = ceil(l1/dz);
z = (1:Nz)*dz;  % space array over EDFA length
Z = (1:N)*dz;   % space array over cavity length
%% array initialization


% Pr = EDFA_ala_p*Pump*EDFA_lam_p/(EDFA_N0*h*c*EDFA_A);           % Pump absorption rate

Pf = zeros(1,N);
Pb = zeros(1,N);
Pp = zeros(1,N);
Ps = zeros(p,Nz);
bp = zeros(1,N);
bs = zeros(p,Nz);
P0 = EDFA_P0*ones(1,N);
N1 = ones(1,N); 
tempf = zeros(1,N);
tempb = zeros(1,N);
temps = zeros(p,Nz);
Ni = -ones(1,Nt);

r13 = (EDFA_ala_p./EDFA_N0)./(h*c./EDFA_lam_p)/(pi*EDFA_wp^2);      % pump absorption rate
r12 = (EDFA_ala_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_ws.^2);      % signal absorption rate
r21 = (EDFA_ale_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_ws.^2);      % stimulated emission rate
A21 = 1/t_f;                                                    % spontaneous emission rate

bs(:,N+nz) = -iso_l/dz;                    % Isolator position
bs(:,N+2*nz) = -cou_l/dz;                  % 3-dB coupler position
bs(:,N+3*nz-10:N+3*nz+14) = -ones(p,25)*aom_l/dz;                  % AOM position
bs(:,N+4*nz) = -wdm_l/dz;                  % WDM position
% ii = find(abs(z-(EDFA_L+3*l1))<dz*1e-10);

Pf(1) = Pumpf;
% PS = zeros(p,Nt);                                               % output power from 3-dB coupler
Pb(1) = Pumpb;
Power = zeros(q,Nt);
Spectrum = zeros(q,p);
%% repetitions

tol = 0.001;
prev_value = 1000;
status = "Not Done";

for r = 1:q
%% Code for Low-Q duration
if status == "Not Done"
    
    l=1;tic
    disp(['OFF-',num2str(r)]);
    t = (1:Nt)*dt;  % time array
    PS = zeros(p,Nt);
    % figure(1)
    for n=1:Nt
    %     tic
        tempf(1) = Pf(1);
        tempf(2:N) = (1+v*dt.*bp(2:N)).*Pf(2:N)-R*(Pf(2:N)-Pf(1:N-1));
        tempb(1) = Pb(1);
        tempb(2:N) = (1 - R + v*dt*fliplr(bp(1:N-1))).*Pb(2:N) + R*Pb(1:N-1);
    %     toc
        temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
        temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
    %     toc
        Pf = tempf;
        Pb = tempb;
        Pp = Pf + fliplr(Pb);
        Ps = temps;
        
        N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(1.-N1));
        bp = -EDFA_ala_p*N1;
        bs(:,1:N) = EDFA_ale_s*(1-N1) - EDFA_ala_s*N1;
        Ps(:,1:N) = Ps(:,1:N) + (EDFA_ale_s*(1-N1)).* P0 * dz;
    %     toc44
        Pf(1) = Pumpf;
        Pb(1) = Pumpb;
        Ni(n) = 1 - 2*N1(N);
        PS(:,n) = Ps(:,N+2*nz-1)-Ps(:,N+2*nz+1);
        temp = sum(PS(:,n));
    
        if abs(prev_value-temp)/prev_value < tol
            status = "Done";
        end
    
        prev_value = temp;
    
    %     disp(n*dt);
        
        if (n*dt > l*1e-7)
            toc
            l=l+1;
    %         subplot(221);semilogy(EDFA_lam_s,PS(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
            subplot(221);plot(Z,Pp)
            subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
            subplot(223);plot(Z,[N1;(1-N1)]);%ylabel('N1/N0');xlabel('length(m)');
            subplot(224);plot(t,sum(PS));%,title('Power'),xlabel('Time'),ylabel('Power');
            pause(1e-3);
            tic
        end
    %     toc
    end
    
    Z_out = Z;
    z_out = z;
    Pp_out = Pp;
    Ps_out = sum(Ps);
    Pout_out = sum(PS);
    inv_out = [N1;(1-N1)];
    t_out = t;
    lambda_s_out = EDFA_lam_s;
    
    subplot(211),plot(t,sum(PS))
    subplot(212),plot(EDFA_lam_s,PS(:,end))
    % figure(2)%,plot((1:Nt)*dt,Ni,'r'),title('Population inversion at end of EDFA'),xlabel('Time'),ylabel('Inversion');
    toc
    % %% Code for High-Q duration
    % l=1;
    % t = (1:Nth)*dt;  % time array
    % disp(['ON-',num2str(r)]);
    % % Ps = zeros(1,Nz);
    % PS = zeros(p,Nth);
    % for n=1:Nth
    %     
    %     tempf(1) = Pf(1);
    %     tempf(2:N) = (1+v*dt.*bp(2:N)).*Pf(2:N)-R*(Pf(2:N)-Pf(1:N-1));
    %     tempb(1) = Pb(1);
    %     tempb(2:N) = (1 - R + v*dt*fliplr(bp(1:N-1))).*Pb(2:N) + R*Pb(1:N-1);
    %     
    %     temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
    %     temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
    % 
    %     Pf = tempf;
    %     Pb = tempb;
    %     Pp = Pf + fliplr(Pb);
    %     Ps = temps;
    %     
    %     N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(1.-N1));
    %     bp = -EDFA_ala_p*N1;
    %     bs(:,1:N) = EDFA_ale_s*(1-N1) - EDFA_ala_s*N1;
    %     Ps(:,1:N) = Ps(:,1:N) + (EDFA_ale_s*(1-N1)).* P0 * dz;
    % %     if(n*dt < tsw)
    % %         bs(:,N+3*nz-10:N+3*nz+14) = -ones(p,25)*AOM(n*dt,2,0.08,tsw)/(dz*4.343);
    % %     end
    % %     if(n*dt>Th-tsw)
    % %         bs(:,N+3*nz-10:N+3*nz+14) = -ones(p,25)*AOM(n*dt-(Th-tsw),0.08,2,tsw)/(dz*4.343);
    % %     end
    %     
    %     Pf(1) = Pumpf;
    %     Pb(1) = Pumpb;
    %     Ni(n) = 1 - 2*N1(N);
    %     PS(:,n) = Ps(:,N+2*nz-1)-Ps(:,N+2*nz);
    % 
    % %     disp(n*dt);
    % %     if (n*dt > l*1e-9)
    % %         l=l+1;
    %         subplot(221);plot(EDFA_lam_s,PS(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
    %         subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
    %         subplot(223);plot(Z,1-2*N1);%ylabel('N1/N0');xlabel('length(m)');
    %         subplot(224);plot(t,sum(PS));%,title('Power'),xlabel('Time'),ylabel('Power');
    %         pause(5e-3)
    % %     end
    % end
    % [a,b]=max(max(PS));toc
    % % figure(1),plot(t,sum(PS)),hold on%title('Power'),xlabel('Time'),ylabel('Power');
    Power(r,:) = sum(PS);
    Spectrum(r,:) = PS(:,end)';
end
end
lam_peak = EDFA_lam_s(Spectrum(r,:)==max(Spectrum(r,:)));
g0 = max(bs((EDFA_lam_s==lam_peak),:));
%Tr = 2*pi/sqrt(v*(g0-Loss)/t_f);
figure()
subplot(211),plot((1:Nt*q)*dt,reshape(Power',1,Nt*q)),legend(int2str((1:q)')),xlabel('Time(s)'),ylabel('Power(W)'),grid on;
subplot(212),plot(EDFA_lam_s*1e9,Spectrum),legend(int2str((1:q)')),xlabel('Wavelength(nm)'),ylabel('Power(W)'),grid on;