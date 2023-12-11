%close all
%clearvars
%clc
% addpath('D:/PhD/M-codes/Functions')

function TDF_CWLaser_0(TDFA_L,TDFA_N0,TDFA_NAs,TDFA_NAp,TDFA_neff,TDFA_rs,TDFA_rp,alpha,lambda_p,lambda_s_start,lambda_s_end,lambda_s_step,tau21,tau32,K1,L,c_l_in,wdm_p_in,iso_l_in,wdm_l_in,Pump_f,Pump_b,dz)
%% initialization of physical constants

nf=1.5; % refractive index of EDFA
c=3e8;  % speed of light
v=c/nf; % velocity of light in EDFA
h=6.625e-34; % Planck's Constant

%% TDF parameters

fiber.L = TDFA_L; %2;                    % Length of gain fiber
fiber.N0 = TDFA_N0; %2.74e26;             % Total dopant concentration in /m^3
fiber.NAs = TDFA_NAs; %0.15;               % Numerical Aperture of the core
fiber.NAp = TDFA_NAp; %0.46;               % Numerical Aperture of the inner clad
fiber.neff = TDFA_neff; %1.5;               % refractive index
fiber.rs = TDFA_rs; %5e-6;                % radius of core in meters
fiber.rp = TDFA_rp; %65e-6;               % radius of clad in meters

fiber.alpha = alpha; %0;

%% Pump and signal channels

fiber.Lp = lambda_p*1e-9; %790e-9;

l1 = lambda_s_start; %1800e-9;
l2 = lambda_s_end; %2100e-9;
dl = lambda_s_step; %2e-9;
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

%% Cavity parameters

L = L; %10;                                        % length of the cavity
l1 = (L-fiber.L)/6;                              % length of the fibre section between 2 components
tr = L/v;                                       % cavity roundtrip time

l_on = 0;
l_off = 56;
al = 40;

c_l = c_l_in; %0;

wdm_p = wdm_p_in/4.343; %1.55/4.343;
iso_l = iso_l_in/4.343; %0.9/4.343;                             % insertion loss of isolator
cou_l = (c_l+3.0)/4.343;
c_out = (c_l+3.0)/4.343;                                % insertion loss of 3-dB coupler
aom_l = l_on/(al*4.343);                                % insertion loss of AOM
tsw = 100e-9;                                    % AOM switching time
wdm_l = wdm_l_in/4.343; %1.1/4.343;                             % insertion loss of WDM

fil_l = 0;
lam_c = 2000e-9;
dL = Inf;
if(dL==Inf)
    tbpf = ones(length(fiber.Ls),1);
else
    sig = dL/(2*sqrt(2*log(2)));
    tbpf = (10^(-fil_l/10))*exp(-0.5*((fiber.Ls - lam_c)/sig).^2);
end
TBPF = -10*log10(tbpf);
TBPF(TBPF>50) = 50;

TBPF = TBPF/4.343;

% r_out = exp(-c_out+TBPF)';
r_out = exp(-c_out);
%% System Parameters

% Pumpf = 1e-3*input('Pump power (mW) = ')*exp(-wdm_p);                                   % Pump power
Pump = Pump_f; %1900;
Pumpf = 1e-3*Pump*exp(-wdm_p);                                   % Pump power
Pumpb = Pump_b; %0;

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

dz = dz; %0.01;
N = fiber.L/dz;

dt = R * dz/v;  % size of time steps
Nt = ceil(T/dt);      % No. of time steps
% Ntl = ceil(Tl/dt);
% Nth = ceil(Th/dt);
t = (1:Nt)*dt;

Nz = ceil(L/dz);      % No. of space steps over the cavity
nz = ceil(l1/dz);
z = (1:Nz)*dz;  % space array over cavity length
Z = (1:N)*dz;   % space array over TDFA length

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
N1 = fiber.N0*ones(1,N);
N2 = zeros(1,N);
N3 = fiber.N0 - N1 - N2;
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
    
    %disp([" TDFA_L = ",TDFA_L]);
    %disp(["TDFA_N0 = ",TDFA_N0]);
    %disp(["TDFA_NAs = ",TDFA_NAs]);
    %disp(["TDFA_NAp = ",TDFA_NAp]);
    %disp(["TDFA_neff = ",TDFA_neff]);
    %disp(["TDFA_rs = ",TDFA_rs]);
    %disp(["TDFA_rp = ",TDFA_rp])
    %disp(["alpha = ",alpha])
    %disp(["lambda_p = ",lambda_p])
    %disp(["lambda_s_start = ",lambda_s_start])
    %disp(["lambda_s_end = ",lambda_s_start])
    %disp(["lambda_s_step = ",lambda_s_start])
    %disp(["tau21 = ",tau21])
    %disp(["tau32 = ",tau21])
    %disp(["K1 = ",K1])
    %disp(["Pump_f = ",Pump_f])

    [Pp_out,Z_out,Spectrum_out,Ps_out,lambda_s,Excited_state_pop_out,Ps_out,N1_out,N2_out,bp_out,bs_out] = TDFA_DCF_0(TDFA_L,TDFA_N0,TDFA_NAs,TDFA_NAp,TDFA_neff,TDFA_rs,TDFA_rp,alpha,lambda_p,lambda_s_start,lambda_s_end,lambda_s_step,tau21,tau32,K1,Pump_f,0,"n",dz);

    pf = Pp_out(1:end-1); %dlmread('Pp.csv');
    ps = Ps_out(1:end-1,1:end-1); %dlmread('Ps.csv');
    n1 = N1_out(1:end-1); %dlmread('N1.csv');
    n2 = N2_out(1:end-1); %dlmread('N2.csv');
    Bp = bp_out(1:end-1); %dlmread('bp.csv');
    Bs = bs_out(1:end-1,1:end-1); %dlmread('bs.csv');
    lams = lambda_s(1:end-1); %dlmread('lams.csv');
    nn = length(n1);
    zz = linspace(fiber.L/nn,fiber.L,nn);
    
    %disp(["size of zz = ",size(zz)]);
    %disp(["size of lams = ",size(lams)]);
    %disp(["size of ps = ",size(ps)]);
    %disp(["size of Z = ",size(Z)]);
    %disp(["size of fiber.Ls(1:end-1) = ",size(fiber.Ls(1:end-1))]);
    %disp(["size of pf = ",size(pf)])
    %disp(["size of Z = ",size(Z)])



    Pf = interp1(zz,pf,Z);
    N1 = interp1(zz,n1,Z);
    N2 = interp1(zz,n2,Z);
    bp = interp1(zz,Bp,Z);

    %disp(["size of Ps(1:end-1,1:N) = ",size(Ps(1:end-1,1:N))])
    %disp(["size of zz = ",size(zz)])
    %disp(["size of lams = ",size(lams)])
    %disp(["size of ps = ",size(ps)])
    %disp(["size of Z = ",size(Z)])
    %disp(["size of fiber.Ls(1:end-1) = ",size(fiber.Ls(1:end-1))])
    %disp(["N = ",N])
    %disp(["size of Bs = ",size(Bs)])

    Ps(1:end-1,1:N) = interp2(zz,lams,ps,Z,fiber.Ls(1:end-1));
    bs(1:end-1,1:N) = interp2(zz,lams,Bs,Z,fiber.Ls(1:end-1));
    Ps(end,1:N) = interp1(zz,ps(end,:),Z);
    bs(end,1:N) = interp1(zz,Bs(end,:),Z);
    N3 = fiber.N0 - N1 - N2;
    
    for r=N+1:Nz
        Ps(:,r) = Ps(:,r-1).*exp(bs(:,r)*dz);
    end
    clear pf ps n1 Bp Bs lams r
elseif(ch=='n')
    Pf(1) = Pumpf;
end

%% Cycles
t0 = 0;
figure()
for r = 1:q
    l=1;tic;
    disp(['ON-',num2str(r)]);
    for n = 1:Nt
        
        tempf(1) = Pf(1);
        tempf(2:N) = (1+v*dt.*bp(2:N)).*Pf(2:N)-R*(Pf(2:N)-Pf(1:N-1));
        
        
        temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
        temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
        
        
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
        PS(:,n) = Ps(:,N+3*nz-1).*r_out';
        
        
        if (n*dt > l*1e-8)
            l=l+1;
            %subplot(221);semilogy(fiber.Ls,PS(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
            %subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
             %         subplot(222);imagesc(z,EDFA_lam_s*1e9,10*log10(1e3*abs(Ps)))
            %subplot(223);plot(Z,(N2-N1)/fiber.N0);%ylabel('N1/N0');xlabel('length(m)');
            %subplot(224);plot(t,sum(PS));%,title('Power'),xlabel('Time'),ylabel('Power');
            %pause(1e-3)
            fprintf('%1.2f us\n',n*dt*1e6)
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

end