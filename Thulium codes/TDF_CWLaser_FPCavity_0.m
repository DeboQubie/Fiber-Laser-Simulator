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
dl = 2e-9;
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

L = 5;                                        % length of the cavity
l1 = (L-fiber.L)/5;                              % length of the fibre section between 2 components
tr = 2*L/v;                                       % cavity roundtrip time

% l_on = 0;
% l_off = 56;
% al = 40;
% 
% c_l = 0;
% 
% wdm_p = 1.55/4.343;
% iso_l = 0.9/4.343;                             % insertion loss of isolator
% cou_l = (c_l+3.0)/4.343;
% c_out = (c_l+3.0)/4.343;                                % insertion loss of 3-dB coupler
% aom_l = l_on/(al*4.343);                                % insertion loss of AOM
% tsw = 100e-9;                                    % AOM switching time
% wdm_l = 1.1/4.343;                             % insertion loss of WDM

fil_l = 0;
lam_c = 2000e-9;
dL = 1000e-9;
sig = dL/(2*sqrt(2*log(2)));
tbpf = (10^(-fil_l/10))*exp(-0.5*((fiber.Ls - lam_c)/sig).^2);

TBPF = -10*log10(tbpf);
TBPF(TBPF>50) = 50;

TBPF = TBPF/4.343;

R1 = 0.99;
R2 = 0.5;
mpc_p = 0.2/4.343;
mpc_l = 0.35/4.343;
mfa_l = 0.25/4.343;


% r_out = exp(-c_out+TBPF)';
r_out = 1 - R2;
%% System Parameters

% Pumpf = 1e-3*input('Pump power (mW) = ')*exp(-wdm_p);                                   % Pump power
Pump = 1900;
Pumpf = 1e-3*Pump*exp(-mpc_p);                                   % Pump power

Tt = 1e-3;
T = 100e-6;
q = Tt/T;

% Rr = 10e3;
% T = 1/Rr;
% Th = 1e-6;
% Tl = T-Th;
% q = 40;
%% Numerical parameters

R = 0.9;        % Stability factor : Courant parameter
% N = 233;         % No. of space steps over EDFA
% dz = EDFA_L/N;  % size of space steps

dz = 0.01;
N = fiber.L/dz;

dt = R * dz/v;  % size of time steps
Nt = ceil(T/dt);      % No. of time steps
% Ntl = ceil(Tl/dt);
% Nth = ceil(Th/dt);
t = (1:Nt)*dt;

Nz = ceil(L/dz);      % No. of space steps over the cavity
nz = ceil(l1/dz);
z = (1:Nz)*dz;  % space array over cavity length
Z = (1:N)*dz;   % space array over EDFA length

tol = 0.01;
%% array initialization


% Pr = EDFA_ala_p*Pump*EDFA_lam_p/(EDFA_N0*h*c*EDFA_A);           % Pump absorption rate

Pp = zeros(1,N);
Psf = zeros(g,Nz);
Psb = zeros(g,Nz);
bp = zeros(1,N);
bs = zeros(g,Nz);
P0 = fiber.PoA*ones(1,N);
N1 = fiber.N0*ones(1,N);
N2 = zeros(1,N);
N3 = fiber.N0 - N1 - N2;
tempp = zeros(1,N);
tempf = zeros(g,Nz);
tempb = zeros(g,Nz);
% Ni = -ones(1,Nt);

% r13 = (EDFA_ala_p./EDFA_N0)./(h*c./EDFA_lam_p)/(pi*EDFA_w^2);      % pump absorption rate
r13 = (fiber.Lp*fiber.sap)/(pi*(fiber.rs^2)*h*c);
% r12 = (EDFA_ala_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_w^2);      % signal absorption rate
% r21 = (EDFA_ale_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_w^2);      % stimulated emission rate
r21 = (fiber.Ls.*fiber.ses)./(pi.*(fiber.rs.^2)*h*c);         %Signal emission rate 
r12 = (fiber.Ls.*fiber.sas)./(pi.*(fiber.rs.^2)*h*c);         %signal absorption rate

A21 = 1/fiber.tau21;                                                    % spontaneous emission rate
A32 = 1/fiber.tau32;

% bs(:,N+nz) = -iso_l/dz;                    % Isolator position
% bs(:,N+2*nz) = -TBPF/dz;
% bs(:,N+2*nz-(al/2):N+2*nz+(al/2)-1) = -TBPF*ones(1,al)/(al*dz);
% bs(:,N+3*nz) = -cou_l/dz;                  % 3-dB coupler position
% bs(:,N+4*nz-(al/2):N+4*nz+(al/2)-1) = -ones(g,al)*aom_l/dz;                  % AOM position
% bs(:,N+5*nz) = -wdm_l/dz;                  % WDM position
% ii = find(abs(z-(EDFA_L+3*l1))<dz*1e-10);

bs(:,nz) = -mpc_l/dz;

bs(:,4*nz) = -mfa_l/dz;

% PSl = zeros(g,ceil(Nth/10));
PS = zeros(g,Nt);
% Power = zeros(q,Nth + length(PSl));                    % output power from 3-dB coupler
Power = zeros(q,Nt);
Spectrum = zeros(q,g);

%% Initial Conditions

ch = input('Start from open-loop initial conditions? y/n [n]: ','s');
if(isempty(ch))
    ch = 'n';
end
if(ch=='y')
    pf = dlmread('Pp.csv');
    ps = dlmread('Ps.csv');
    n1 = dlmread('N1.csv');
    n2 = dlmread('N2.csv');
    Bp = dlmread('bp.csv');
    Bs = dlmread('bs.csv');
    lams = dlmread('lams.csv');
    nn = length(n1);
    zz = linspace(fiber.L/nn,fiber.L,nn);
    
    Pp = interp1(zz,pf,Z);
    N1 = interp1(zz,n1,Z);
    N2 = interp1(zz,n2,Z);
    bp = interp1(zz,Bp,Z);
    Psf(1:end-1,2*nz+(1:N)) = interp2(zz,lams,ps,Z,fiber.Ls(1:end-1));
    bs(1:end-1,2*nz+(1:N)) = interp2(zz,lams,Bs,Z,fiber.Ls(1:end-1));
    Psf(end,2*nz+(1:N)) = interp1(zz,ps(end,:),Z);
    bs(end,2*nz+(1:N)) = interp1(zz,Bs(end,:),Z);
    N3 = fiber.N0 - N1 - N2;
    
%     for r=N+1:Nz
%         Ps(:,r) = Ps(:,r-1).*exp(bs(:,r)*dz);
%     end
    clear pf ps n1 Bp Bs lams r
elseif(ch=='n')
    Pp(1) = Pumpf;
end

%%
t0 = 0;
figure()
for r = 1:q
    l=1;tic;
    disp(['ON-',num2str(r)]);
    for n = 1:Nt
        
        tempp(1) = Pp(1);
        tempp(2:N) = (1+v*dt.*bp(2:N)).*Pp(2:N)-R*(Pp(2:N)-Pp(1:N-1));
        
        tempf(:,1) = R1*Psb(:,Nz);
        tempf(:,2:Nz) = (1+v*dt.*bs(:,2:Nz)).*Psf(:,2:Nz)-R*(Psf(:,2:Nz)-Psf(:,1:Nz-1));
        tempb(:,1) = R2*Psf(:,Nz);
        tempb(:,2:Nz) = (1+v*dt.*fliplr(bs(:,2:Nz))).*Psb(:,2:Nz)-R*(Psb(:,2:Nz)-Psb(:,1:Nz-1));
        
        
%         temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
%         temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
%         
        
        Pp = tempp;
        Psf = tempf;
        Psb = tempb;
        Ps = Psf + fliplr(Psb);
%         Ps = temps;
        
        CR = fiber.K1*N1.*N3 - fiber.K2*N2.*N2;
        
        N3 = fiber.N0 - N1 - N2;
        N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,2*nz+(1:N)))).*N1 + (A21+(r21'*Ps(:,2*nz+(1:N)))).*(N2) - CR);
        N2 = N2 + dt*(-(A21 + (r21'*Ps(:,2*nz+(1:N)))).*N2 + (r12'*Ps(:,2*nz+(1:N))).*N1 + A32*N3 + 2*CR);
        %     N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(1.-N1));
        bp = -fiber.sap*N1;
        bs(:,2*nz+(1:N)) = fiber.ses*(N2) - fiber.sas*N1;
        Psf(:,2*nz+(1:N)) = Psf(:,2*nz+(1:N)) + (fiber.ses*(N2)).* P0 * dz;
        Psb(:,2*nz+(1:N)) = Psb(:,2*nz+(1:N)) + (fiber.ses*(fliplr(N2))).* P0 * dz;
        
        
        Pp(1) = Pumpf;
        %     Ni(n) = 1 - 2*N1(N);
        PS(:,n) = Psf(:,Nz).*r_out';
        
        
        if (n*dt > l*1e-9)
            l=l+1;
            subplot(221);semilogy(fiber.Ls,PS(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
            subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
            %         subplot(222);imagesc(z,EDFA_lam_s*1e9,10*log10(1e3*abs(Ps)))
            subplot(223);plot(Z,(N2-N1)/fiber.N0);%ylabel('N1/N0');xlabel('length(m)');
            subplot(224);plot(t,sum(PS));%,title('Power'),xlabel('Time'),ylabel('Power');
            pause(1e-3)
%             fprintf('%1.2f us\n',n*dt*1e6)
        end
    end
%     [a,b]=max(max(PS));
    toc
    Spectrum(r,:) = PS(:,end)';
    Power(r,:) = sum(PS);
    plot(t0 + t,Power(r,:),'b','linewidth',3),pause(1e-3)
    hold on
    t0 = r*T;
end
figure()
semilogy(fiber.Ls,Spectrum,'linewidth',3)
G = trapz(Z,bs(:,1:N)');
figure()
plot(fiber.Ls,4.343*G,'linewidth',3)
figure()
plot(Z,[N1;N2;N3]/fiber.N0,'linewidth',3)
