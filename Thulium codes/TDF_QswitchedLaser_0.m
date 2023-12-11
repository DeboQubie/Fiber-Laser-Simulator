close all
clearvars
clc
% addpath('D:/PhD/M-codes/Functions')
%% initialization of physical constants

nf=1.5; % refractive index of EDFA
c=3e8;  % speed of light
v=c/nf; % velocity of light in EDFA
h=6.625e-34; % Planck's Constant

%% TDF parameters

fiber.L = 2;                    % Length of gain fiber
fiber.N0 = 2.74e26;             % Total dopant concentration in /m^3
fiber.NAs = 0.15;               % Numerical Aperture of the core
fiber.NAp = 0.46;               % Numerical Aperture of the inner clad
fiber.neff = 1.5;               % refractive index
fiber.rs = 5e-6;                % radius of core in meters
fiber.rp = 65e-6;               % radius of clad in meters

fiber.alpha = 0;

%% Pump and signal channels

fiber.Lp = 790e-9;

l1 = 1800e-9;
l2 = 2100e-9;
dl = 1e-9;
fiber.Ls = (l1:dl:l2)';

g = length(fiber.Ls);

[fiber.wp,~] = mfd(fiber.Lp,fiber.rp,fiber.NAp);
[fiber.ws,fiber.gammas] = mfd(fiber.Ls,fiber.rs,fiber.NAs);
fiber.gammap = (fiber.rs/fiber.rp)^2;

%% Spectroscopic parameters, absorption and emission cross sections

fiber.tau21 = 334.7e-6;         % Lifetime of level 2 in sec
fiber.tau32 = 14.2e-6;          % Lifetime of level 3 in sec

fiber.K1 = 1.25e-22;            % K3212 cross relaxation coefficient
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

%% Cavity parameters

L = 10;                                        % length of the cavity
l1 = (L-fiber.L)/6;                              % length of the fibre section between 2 components
tr = L/v;                                       % cavity roundtrip time

l_on = 3;
l_off = 56;
al = 40;

c_l = 0;

wdm_p = 1.55/4.343;
iso_l = 0.9/4.343;                             % insertion loss of isolator
cou_l = (c_l+3.35)/4.343;
c_out = (c_l+3.35)/4.343;                                % insertion loss of 3-dB coupler
aom_l = l_off/(al*4.343);                                % insertion loss of AOM
tsw = 100e-9;                                    % AOM switching time
wdm_l = 1.6/4.343;                             % insertion loss of WDM

fil_l = 2;
lam_c = 2000e-9;
dL = 1000e-9;
sig = dL/(2*sqrt(2*log(2)));
tbpf = (10^(-fil_l/10))*exp(-0.5*((fiber.Ls - lam_c)/sig).^2);

TBPF = -10*log10(tbpf);
TBPF(TBPF>50) = 50;

TBPF = TBPF/4.343;

% r_out = exp(-c_out+TBPF)';
r_out = exp(-c_out);
%% System Parameters

% Pumpf = 1e-3*input('Pump power (mW) = ')*exp(-wdm_p);                                   % Pump power
Pump = 1584;
Pumpf = 1e-3*Pump*exp(-wdm_p);                                   % Pump power
Pumpb = 0;

Rr = 10e3;
T = 1/Rr;
Th = 1e-6;
Tl = T-Th;
q = 40;
%% Numerical parameters

R = 0.9;        % Stability factor : Courant parameter
% N = 233;         % No. of space steps over EDFA
% dz = EDFA_L/N;  % size of space steps

dz = 0.02;
N = fiber.L/dz;

dt = R * dz/v;  % size of time steps
Nt = ceil(T/dt);      % No. of time steps
Ntl = ceil(Tl/dt);
Nth = ceil(Th/dt);

Nz = ceil(L/dz);      % No. of space steps over the cavity
nz = ceil(l1/dz);
z = (1:Nz)*dz;  % space array over cavity length
Z = (1:N)*dz;   % space array over EDFA length

tol = 0.01;
%% array initialization


% Pr = EDFA_ala_p*Pump*EDFA_lam_p/(EDFA_N0*h*c*EDFA_A);           % Pump absorption rate

Pf = zeros(1,N);
Pb = zeros(1,N);
Pp = zeros(1,N);
Ps = zeros(g,Nz);
bp = zeros(1,N);
bs = zeros(g,Nz);
P0 = fiber.PoA*ones(1,N);
N1 = ones(1,N); 
tempf = zeros(1,N);
tempb = zeros(1,N);
temps = zeros(g,Nz);
% Ni = -ones(1,Nt);

% r13 = (EDFA_ala_p./EDFA_N0)./(h*c./EDFA_lam_p)/(pi*EDFA_w^2);      % pump absorption rate
r13 = (fiber.Lp*fiber.sap)/(pi*(fiber.rs^2)*h*c);
% r12 = (EDFA_ala_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_w^2);      % signal absorption rate
% r21 = (EDFA_ale_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_w^2);      % stimulated emission rate
r21 = (fiber.Ls.*fiber.ses)./(pi.*(fiber.rs.^2)*h*c);         %Signal emission rate 
r12 = (fiber.Ls.*fiber.sas)./(pi.*(fiber.rs.^2)*h*c);         %signal absorption rate

A21 = 1/fiber.tau21;                                                    % spontaneous emission rate
A32 = 1/fiber.tau32;

bs(:,N+nz) = -iso_l/dz;                    % Isolator position
% bs(:,N+2*nz) = -TBPF/dz;
bs(:,N+2*nz-(al/2):N+2*nz+(al/2)-1) = -TBPF*ones(1,al)/(al*dz);
bs(:,N+3*nz) = -cou_l/dz;                  % 3-dB coupler position
bs(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) = -ones(g,al)*aom_l/dz;                  % AOM position
bs(:,N+5*nz) = -wdm_l/dz;                  % WDM position
% ii = find(abs(z-(EDFA_L+3*l1))<dz*1e-10);

PSl = zeros(g,ceil(Nth/10));
PSh = zeros(g,Nth);
Power = zeros(q,Nth + length(PSl));                    % output power from 3-dB coupler
% Power = zeros(q,Nth);
Spectrum = zeros(q,g);
%% Initial Conditions

initial_condition = "Closed";

%if initial_condition == "Open"
    
%    TDFA_DCF_0(fiber.L,fiber.N0,fiber.NAs,fiber.NAp,fiber.neff,fiber.rs,fiber.rp,fiber.alpha,fiber.Lp*1e9,fiber.tau21,fiber.tau32,fiber.K1,Pump, 0,"y");

%else

%    TDFRL(fiber.L,fiber.N0,fiber.rs,fiber.NAs,fiber.alpha,Pump,fiber.Lp*1e9,fiber.tau21,fiber.tau32,"y");

%end

pf = dlmread('Pp.csv');
ps = dlmread('Ps.csv');
ps = ps';
n1 = dlmread('N1.csv');
n2 = dlmread('N2.csv');
Bp = dlmread('bp.csv');
Bs = dlmread('bs.csv');
lams = dlmread('lams.csv');
nn = length(n1);
zz = linspace(fiber.L/nn,fiber.L,nn);


disp(size(zz));
disp(size(lams));
disp(size(ps));
disp(size(Z'));
disp(size(fiber.Ls(1:end-1)));

Pf = interp1(zz,pf,Z);
N1 = interp1(zz,n1,Z);
N2 = interp1(zz,n2,Z);
bp = interp1(zz,Bp,Z);
Ps(1:end-1,1:N) = interp2(zz,lams,ps,Z,fiber.Ls(1:end-1));
bs(1:end-1,1:N) = interp2(zz,lams,Bs,Z,fiber.Ls(1:end-1));
Ps(end,1:N) = interp1(zz,ps(end,:),Z);
bs(end,1:N) = interp1(zz,Bs(end,:),Z);
N3 = fiber.N0 - N1 - N2;

for r=N+1:Nz
    Ps(:,r) = Ps(:,r-1).*exp(bs(:,r)*dz);
end
clear pf ps n1 Bp Bs lams r
%% repetitions

for r = 1:q
%% Code for High-Q duration
l=1;tic;
t = (1:Nth)*dt;  % time array
disp(['ON-',num2str(r)]);
% Ps = zeros(1,Nz);
PSh = zeros(g,Nth);
for n=1:Nth
    
    tempf(1) = Pf(1);
    tempf(2:N) = (1+v*dt.*bp(2:N)).*Pf(2:N)-R*(Pf(2:N)-Pf(1:N-1));
%     tempb(1) = Pb(1);
%     tempb(2:N) = (1 - R + v*dt*fliplr(bp(1:N-1))).*Pb(2:N) + R*Pb(1:N-1);

    temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
    temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
    
%     temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
%     temps(:,2:N) = (1+v*dt.*bs(:,2:N) - R).*Ps(:,2:N) + R*Ps(:,1:N-1);
%     
%     temps(:,N+1:N+nz-1) = Ps(:,N:N+nz-2);
%     temps(:,N+nz) = Ps(:,N+nz-1).*exp(bs(:,N+nz)*dz);
%     temps(:,N+nz+1:N+2*nz-1) = Ps(:,N+nz:N+2*nz-2);
%     temps(:,N+2*nz) = Ps(:,N+2*nz-1).*exp(bs(:,N+2*nz)*dz);
%     temps(:,N+2*nz+1:N+3*nz-1) = Ps(:,N+2*nz:N+3*nz-2);
%     temps(:,N+3*nz) = Ps(:,N+3*nz-1).*exp(bs(:,N+3*nz)*dz);
%     temps(:,N+3*nz+1:N+4*nz-(al/2)-1) = Ps(:,N+3*nz:N+4*nz-(al/2)-2);
%     temps(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) = (1+v*dt.*bs(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) - R).*Ps(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) + R*Ps(:,N+4*nz-(al/2)-1:N+4*nz+(al/2)-2);
%     temps(:,N+4*nz+(al/2):N+5*nz-1) = Ps(:,N+4*nz+(al/2)-1:N+5*nz-2);
%     temps(:,N+5*nz) = Ps(:,N+5*nz-1).*exp(bs(:,N+5*nz)*dz);
%     temps(:,N+5*nz+1:Nz) = Ps(:,N+5*nz:Nz-1);

    Pf = tempf;
    Pb = tempb;
    Pp = Pf + fliplr(Pb);
    Ps = temps;
    
    CR = fiber.K1*N1.*N3 - fiber.K2*N2.*N2;
    
    N3 = fiber.N0 - N1 - N2;
    N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(N2) - CR);
    N2 = N2 + dt*(-(A21 + (r21'*Ps(:,1:N))).*N2 + (r12'*Ps(:,1:N)).*N1 + A32*N3 + 2*CR);
%     N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(1.-N1));
    bp = -fiber.sap*N1;
    bs(:,1:N) = fiber.ses*(N2) - fiber.sas*N1;
    Ps(:,1:N) = Ps(:,1:N) + (fiber.ses*(N2)).* P0 * dz;
    if(n*dt < tsw)
        bs(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) = ones(g,al)*AOMw(n*dt,(l_off/al),(l_on/al),tsw)/(dz*4.343);
    end
    if(n*dt>Th-tsw)
        bs(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) = ones(g,al)*AOMw(n*dt-(Th-tsw),(l_on/al),(l_off/al),tsw)/(dz*4.343);
    end
    
    Pf(1) = Pumpf;
    Pb(1) = Pumpb;
%     Ni(n) = 1 - 2*N1(N);
    PSh(:,n) = Ps(:,N+3*nz-1).*r_out';

%     disp(n*dt);
    if (n*dt > l*1e-9)
        l=l+1;
        subplot(221);semilogy(fiber.Ls,PSh(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
        subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
%         subplot(222);imagesc(z,EDFA_lam_s*1e9,10*log10(1e3*abs(Ps)))        
        subplot(223);plot(Z,(N2-N1)/fiber.N0);%ylabel('N1/N0');xlabel('length(m)');
        subplot(224);plot(t,sum(PSh));%,title('Power'),xlabel('Time'),ylabel('Power');
        pause(1e-3)
    end
end
[a,b]=max(max(PSh));toc
% figure(1),plot(t,sum(PS)),hold on%title('Power'),xlabel('Time'),ylabel('Power');
% Power(r,:) = sum(PS);
Spectrum(r,:) = PSh(:,b)';
disp(['Peak power = ',num2str(max(sum(PSh))),' W'])
if(r>2)
    err = abs(Power(r-2,1:Nth)-sum(PSh))./Power(r-2,1:Nth);
end
if(r>2 && max(err)<tol)
    break
end
%% Code for Low-Q duration
l=1;
disp(['OFF-',num2str(r)]);
t = (1:Ntl)*dt;  % time array
PSl = zeros(g,ceil(Nth/10));
% figure(1)
for n=1:Ntl
    
    tempf(1) = Pf(1);
    tempf(2:N) = (1+v*dt.*bp(2:N)).*Pf(2:N)-R*(Pf(2:N)-Pf(1:N-1));
%     tempb(1) = Pb(1);
%     tempb(2:N) = (1 - R + v*dt*fliplr(bp(1:N-1))).*Pb(2:N) + R*Pb(1:N-1);

    temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
    temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
    
%     temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
%     temps(:,2:N) = (1+v*dt.*bs(:,2:N) - R).*Ps(:,2:N) + R*Ps(:,1:N-1);
%     
%     temps(:,N+1:N+nz-1) = Ps(:,N:N+nz-2);
%     temps(:,N+nz) = Ps(:,N+nz-1).*exp(bs(:,N+nz)*dz);
%     temps(:,N+nz+1:N+2*nz-1) = Ps(:,N+nz:N+2*nz-2);
%     temps(:,N+2*nz) = Ps(:,N+2*nz-1).*exp(bs(:,N+2*nz)*dz);
%     temps(:,N+2*nz+1:N+3*nz-1) = Ps(:,N+2*nz:N+3*nz-2);
%     temps(:,N+3*nz) = Ps(:,N+3*nz-1).*exp(bs(:,N+3*nz)*dz);
%     temps(:,N+3*nz+1:N+4*nz-(al/2)-1) = Ps(:,N+3*nz:N+4*nz-(al/2)-2);
%     temps(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) = (1+v*dt.*bs(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) - R).*Ps(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) + R*Ps(:,N+4*nz-(al/2)-1:N+4*nz+(al/2)-2);
%     temps(:,N+4*nz+(al/2):N+5*nz-1) = Ps(:,N+4*nz+(al/2)-1:N+5*nz-2);
%     temps(:,N+5*nz) = Ps(:,N+5*nz-1).*exp(bs(:,N+5*nz)*dz);
%     temps(:,N+5*nz+1:Nz) = Ps(:,N+5*nz:Nz-1);
    
    Pf = tempf;
    Pb = tempb;
    Pp = Pf + fliplr(Pb);
    Ps = temps;
    
    CR = fiber.K1*N1.*N3 - fiber.K2*N2.*N2;
    
    N3 = fiber.N0 - N1 - N2;
    N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(N2) - CR);
    N2 = N2 + dt*(-(A21 + (r21'*Ps(:,1:N))).*N2 + (r12'*Ps(:,1:N)).*N1 + A32*N3 + 2*CR);
%     N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(1.-N1));
    bp = -fiber.sap*N1;
    bs(:,1:N) = fiber.ses*(N2) - fiber.sas*N1;
    Ps(:,1:N) = Ps(:,1:N) + (fiber.ses*(N2)).* P0 * dz;
    
    Pf(1) = Pumpf;
    Pb(1) = Pumpb;
%     Ni(n) = 1 - 2*N1(N);
    if(n<=Nth/10)
        PSl(:,n) = Ps(:,N+3*nz-1).*r_out';
%     else if(r==q)
%             break
%         end
    end
    
%     disp(n*dt);
    
%     if (n*dt > l*1e-9)
%         l=l+1;
%         subplot(221);semilogy(fiber.Ls,PSl(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
% %         subplot(221);plot(Z,Pp);
%         subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
%         subplot(223);plot(Z,(N2-N1)/fiber.N0);%ylabel('N1/N0');xlabel('length(m)');
%         subplot(224);plot(t,sum(PSl));%,title('Power'),xlabel('Time'),ylabel('Power');
%         pause(1e-3);
%     end
end
toc
% figure(2)%,plot((1:Nt)*dt,Ni,'r'),title('Population inversion at end of EDFA'),xlabel('Time'),ylabel('Inversion');
Power(r,:) = [sum(PSh) sum(PSl)];
end
% Power(r+1,:) = [sum(PSh) zeros(1,Nth/10)];

% dlmwrite('Ps.csv',Ps)
% dlmwrite('N1.csv',N1)
% dlmwrite('Pf.csv',Pf)
% dlmwrite('bp.csv',bp)
% dlmwrite('bs.csv',bs)
% dlmwrite('Power.csv',Power)
% dlmwrite('Spectrum.csv',Spectrum)
t = (1:(1.1*Nth))*dt;
figure()
subplot(211),plot(t,Power),legend(int2str((1:q)')),xlabel('Time (s)'),ylabel('Power (W)'),title('Pulse Output'),grid on;
subplot(212),plot(EDFA_lam_s*1e9,10*log10(Spectrum*1e3)),legend(int2str((1:q)')),xlabel('Wavelength (nm)'),ylabel('Power (dBm)'),title('Spectrum of Peak'),grid on;