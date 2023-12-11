clear all
close all
clc
%% initialization of physical constants

nf=1.5; % refractive index of EDFA
c=3e8;  % speed of light
v=c/nf; % velocity of light in EDFA
h=6.625e-34; % Planck's Constant
%% EDFA parameters
EDFA_N0 = 2.0e25;                               % Doping concentration
EDFA_L = 3.5;                                    % length of EDFA
EDFA_w = 4.4e-6/2;                               % EDFA Mode field radius
NA = 0.29;

Data = dlmread('HG980.csv');

EDFA_lam_p = Data(1,1)*1e-9;                            % Pump wavelength
[EDFA_wp,gamma_p] = mfd(EDFA_lam_p,EDFA_w,NA);
EDFA_ala_p = Data(1,2)*EDFA_N0*gamma_p*1e-25;                   % abs.coeff at 980

EDFA_lam_s = Data(2:end,1)*1e-9;                  % signal wavelengths 1530nm,1540nm
[EDFA_ws,gamma_s] = mfd(EDFA_lam_s,EDFA_w,NA);

EDFA_ala_s = gamma_s.*Data(2:end,2)*EDFA_N0*1e-25;  % abs. coeffs at 1530,1540
EDFA_ale_s = gamma_s.*Data(2:end,3)*EDFA_N0*1e-25;  % emission coeffs at 1530,1540
EDFA_d_lam = EDFA_lam_s(2)-EDFA_lam_s(1);        % Noise photon bandwidth
EDFA_P0 = h*c^2./EDFA_lam_s.^3*EDFA_d_lam*2;   % spectral emission

t_f = 10e-3;                                     % flourescence lifetime of upper lasing level
p = length(EDFA_lam_s);

EDFA_A = pi*EDFA_w^2;                            % EDFA cross-section area
clear Data
%% Cavity parameters

L = 9.5;                                        % length of the cavity
l1 = (L-EDFA_L)/5;                              % length of the fibre section between 2 components
tr = L/v;                                       % cavity roundtrip time

l_on = 9.5;
l_off = 52.5;
al = 40;

c_l = 2;

wdm_p = 1.55/4.343;
iso_l = 0.9/4.343;                             % insertion loss of isolator
cou_l = (c_l+3.35)/4.343;
c_out = (c_l+3.35)/4.343;                                % insertion loss of 3-dB coupler
aom_l = l_off/(al*4.343);                                % insertion loss of AOM
tsw = 150e-9;                                    % AOM switching time
wdm_l = 1.6/4.343;                             % insertion loss of WDM

r_out = exp(-c_out);
%% System Parameters

Pumpf = 1e-3*input('Pump power (mW) = ')*exp(-wdm_p);                                   % Pump power
Pumpb = 0;

Rr = 10e3;
T = 1/Rr;
Th = 1e-6;
Tl = T-Th;
q = 40;
%% Numerical parameters

R = 0.9;        % Stability factor : Courant parameter
N = 233;         % No. of space steps over EDFA
dz = EDFA_L/N;  % size of space steps
dt = R * dz/v;  % size of time steps
Nt = ceil(T/dt);      % No. of time steps
Ntl = ceil(Tl/dt);
Nth = ceil(Th/dt);

Nz = ceil(L/dz);      % No. of space steps over the cavity
nz = ceil(l1/dz);
z = (1:Nz)*dz;  % space array over EDFA length
Z = (1:N)*dz;   % space array over cavity length

tol = 0.01;
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
% Ni = -ones(1,Nt);

r13 = (EDFA_ala_p./EDFA_N0)./(h*c./EDFA_lam_p)/(pi*EDFA_w^2);      % pump absorption rate
r12 = (EDFA_ala_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_w^2);      % signal absorption rate
r21 = (EDFA_ale_s./EDFA_N0)./(h*c./EDFA_lam_s)./(pi*EDFA_w^2);      % stimulated emission rate
A21 = 1/t_f;                                                    % spontaneous emission rate

bs(:,N+nz) = -iso_l/dz;                    % Isolator position
bs(:,N+2*nz) = -cou_l/dz;                  % 3-dB coupler position
bs(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) = -ones(p,al)*aom_l/dz;                  % AOM position
bs(:,N+4*nz) = -wdm_l/dz;                  % WDM position
ii = find(abs(z-(EDFA_L+3*l1))<dz*1e-10);

PSl = zeros(p,Nth/10);
PSh = zeros(p,Nth);
Power = zeros(q,1.1*Nth);                    % output power from 3-dB coupler
% Power = zeros(q,Nth);
Spectrum = zeros(q,p);
%% Initial Conditions

pf = dlmread('Pp.csv');
ps = dlmread('Ps.csv');
n1 = dlmread('N1.csv');
Bp = dlmread('bp.csv');
Bs = dlmread('bs.csv');
lams = dlmread('lams.csv');
nn = length(n1);
zz = linspace(EDFA_L/nn,EDFA_L,nn);

Pf = interp1(zz,pf,Z);
N1 = interp1(zz,n1,Z);
bp = interp1(zz,Bp,Z);
Ps(1:end-1,1:N) = interp2(zz,lams,ps,Z,EDFA_lam_s(1:end-1));
bs(1:end-1,1:N) = interp2(zz,lams,Bs,Z,EDFA_lam_s(1:end-1));
Ps(end,1:N) = interp1(zz,ps(end,:),Z);
bs(end,1:N) = interp1(zz,Bs(end,:),Z);

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
PSh = zeros(p,Nth);
for n=1:Nth
    
    tempf(1) = Pf(1);
    tempf(2:N) = (1+v*dt.*bp(2:N)).*Pf(2:N)-R*(Pf(2:N)-Pf(1:N-1));
%     tempb(1) = Pb(1);
%     tempb(2:N) = (1 - R + v*dt*fliplr(bp(1:N-1))).*Pb(2:N) + R*Pb(1:N-1);

%     temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
%     temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
    
    temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
    temps(:,2:N) = (1+v*dt.*bs(:,2:N) - R).*Ps(:,2:N) + R*Ps(:,1:N-1);
    
    temps(:,N+1:N+nz-1) = Ps(:,N:N+nz-2);
    temps(:,N+nz) = Ps(:,N+nz-1).*exp(bs(:,N+nz)*dz);
    temps(:,N+nz+1:N+2*nz-1) = Ps(:,N+nz:N+2*nz-2);
    temps(:,N+2*nz) = Ps(:,N+2*nz-1).*exp(bs(:,N+2*nz)*dz);
    temps(:,N+2*nz+1:N+3*nz-(al/2)-1) = Ps(:,N+2*nz:N+3*nz-(al/2)-2);
    temps(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) = (1+v*dt.*bs(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) - R).*Ps(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) + R*Ps(:,N+3*nz-(al/2)-1:N+3*nz+(al/2)-2);
    temps(:,N+3*nz+(al/2):N+4*nz-1) = Ps(:,N+3*nz+(al/2)-1:N+4*nz-2);
    temps(:,N+4*nz) = Ps(:,N+4*nz-1).*exp(bs(:,N+4*nz)*dz);
    temps(:,N+4*nz+1:Nz) = Ps(:,N+4*nz:Nz-1);

    Pf = tempf;
    Pb = tempb;
    Pp = Pf + fliplr(Pb);
    Ps = temps;
    
    N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(1.-N1));
    bp = -EDFA_ala_p*N1;
    bs(:,1:N) = EDFA_ale_s*(1-N1) - EDFA_ala_s*N1;
    Ps(:,1:N) = Ps(:,1:N) + (EDFA_ale_s*(1-N1)).* P0 * dz;
    if(n*dt < tsw)
        bs(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) = -ones(p,al)*AOM(n*dt,(l_off/al),(l_on/al),tsw)/(dz*4.343);
    end
    if(n*dt>Th-tsw)
        bs(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) = -ones(p,al)*AOM(n*dt-(Th-tsw),(l_on/al),(l_off/al),tsw)/(dz*4.343);
    end
    
    Pf(1) = Pumpf;
    Pb(1) = Pumpb;
%     Ni(n) = 1 - 2*N1(N);
    PSh(:,n) = Ps(:,N+2*nz-1)*r_out;

%     disp(n*dt);
    if (n*dt > l*1e-9)
        l=l+1;
        subplot(221);semilogy(EDFA_lam_s,PSh(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
        subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
%         subplot(222);imagesc(z,EDFA_lam_s*1e9,10*log10(1e3*abs(Ps)))
        subplot(223);plot(Z,1-2*N1);%ylabel('N1/N0');xlabel('length(m)');
        subplot(224);plot(t,sum(PSh));%,title('Power'),xlabel('Time'),ylabel('Power');
        pause(5e-3)
    end
end
[a,b]=max(max(PSh));toc
% figure(1),plot(t,sum(PS)),hold on%title('Power'),xlabel('Time'),ylabel('Power');
% Power(r,:) = sum(PS);
Spectrum(r,:) = PSh(:,b)';
disp(['Peak power = ',num2str(max(sum(PSh)))])
if(r>1)
    err = abs(Power(r-1,1:Nth)-sum(PSh))./Power(r-1,1:Nth);
end
if(r>1 && max(err)<tol)
    break
end
%% Code for Low-Q duration
l=1;
disp(['OFF-',num2str(r)]);
t = (1:Ntl)*dt;  % time array
PSl = zeros(p,Nth/10);
% figure(1)
for n=1:Ntl
    
    tempf(1) = Pf(1);
    tempf(2:N) = (1+v*dt.*bp(2:N)).*Pf(2:N)-R*(Pf(2:N)-Pf(1:N-1));
%     tempb(1) = Pb(1);
%     tempb(2:N) = (1 - R + v*dt*fliplr(bp(1:N-1))).*Pb(2:N) + R*Pb(1:N-1);

%     temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
%     temps(:,2:Nz) = (1+v*dt.*bs(:,2:Nz) - R).*Ps(:,2:Nz) + R*Ps(:,1:Nz-1);
    
    temps(:,1) = (1 + v*dt*bs(:,1) - R).*Ps(:,1) + R * Ps(:,Nz);
    temps(:,2:N) = (1+v*dt.*bs(:,2:N) - R).*Ps(:,2:N) + R*Ps(:,1:N-1);
    
    temps(:,N+1:N+nz-1) = Ps(:,N:N+nz-2);
    temps(:,N+nz) = Ps(:,N+nz-1).*exp(bs(:,N+nz)*dz);
    temps(:,N+nz+1:N+2*nz-1) = Ps(:,N+nz:N+2*nz-2);
    temps(:,N+2*nz) = Ps(:,N+2*nz-1).*exp(bs(:,N+2*nz)*dz);
    temps(:,N+2*nz+1:N+3*nz-(al/2)-1) = Ps(:,N+2*nz:N+3*nz-(al/2)-2);
    temps(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) = (1+v*dt.*bs(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) - R).*Ps(:,N+3*nz-(al/2):N+3*nz+(al/2)-1) + R*Ps(:,N+3*nz-(al/2)-1:N+3*nz+(al/2)-2);
    temps(:,N+3*nz+(al/2):N+4*nz-1) = Ps(:,N+3*nz+(al/2)-1:N+4*nz-2);
    temps(:,N+4*nz) = Ps(:,N+4*nz-1).*exp(bs(:,N+4*nz)*dz);
    temps(:,N+4*nz+1:Nz) = Ps(:,N+4*nz:Nz-1);
    
    Pf = tempf;
    Pb = tempb;
    Pp = Pf + fliplr(Pb);
    Ps = temps;
    
    N1 = N1 + dt*(-(r13'*Pp+(r12'*Ps(:,1:N))).*N1 + (A21+(r21'*Ps(:,1:N))).*(1.-N1));
    bp = -EDFA_ala_p*N1;
    bs(:,1:N) = EDFA_ale_s*(1-N1) - EDFA_ala_s*N1;
    Ps(:,1:N) = Ps(:,1:N) + (EDFA_ale_s*(1-N1)).* P0 * dz;
    
    Pf(1) = Pumpf;
    Pb(1) = Pumpb;
%     Ni(n) = 1 - 2*N1(N);
    if(n<=Nth/10)
        PSl(:,n) = Ps(:,N+2*nz-1)*r_out;
%     else if(r==q)
%             break
%         end
    end
    
%     disp(n*dt);
    
%     if (n*dt > l*1e-9)
%         l=l+1;
% %         subplot(221);semilogy(EDFA_lam_s,PSl(:,n)');%ylabel('gain coefficient(dB/m)');xlabel('length(m)');
%         subplot(221);plot(Z,Pp);
%         subplot(222);plot(z,sum(Ps));%ylabel('Signal(W)');xlabel('length(m)');
%         subplot(223);plot(Z,1-2*N1);%ylabel('N1/N0');xlabel('length(m)');
% %         subplot(224);plot(t,sum(PSl));%,title('Power'),xlabel('Time'),ylabel('Power');
%         pause(5e-3);
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
subplot(211),plot(t,Power),legend(int2str((1:q)')),xlabel('Time (s)'),ylabel('Power (W)'),grid on;
subplot(212),plot(EDFA_lam_s*1e9,10*log10(Spectrum*1e3)),legend(int2str((1:q)')),xlabel('Wavelength(nm)'),ylabel('Power(dBm)'),grid on;