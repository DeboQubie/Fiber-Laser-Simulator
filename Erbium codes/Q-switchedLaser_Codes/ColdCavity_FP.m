close all
clearvars
clc
addpath('D:/PhD/M-codes/Functions')
%% initialization of physical constants

nf=1.5; % refractive index of EDFA
c=3e8;  % speed of light
v=c/nf; % velocity of light in EDFA
h=6.625e-34; % Planck's Constant
%% cavity

R1 = 0.99;
R2 = 0.5;

L = 2.5;

tr = 2*L/v;

r1 = sqrt(R1);
r2 = sqrt(R2);
tout = sqrt(1-R2);
tin = sqrt(1-R1);

%% Frequency/wavelength axes

lam_c = 1550e-9;
f_c = c/lam_c;

fsr = 1/tr;

df = fsr/5;
Fspan = 1e9;
f = -Fspan/2:df:Fspan/2;
F = f + f_c;

lam_s = c./F;

g = length(f);

%% System parameters

T = 100e-6;
q = 10;

%% Numerical parameters

R = 0.9;
dz = 0.02;

N = L/dz;

dt = R*dz/v;
Nt = ceil(T/dt);

z = (1:N)*dz;
t = (1:Nt)*dt;

%% Array initialization

Es = zeros(g,N);
Ef = zeros(g,N);
Eb = zeros(g,N);
bs = zeros(g,N);

tempf = zeros(g,N);
tempb = zeros(g,N);


Field_f = zeros(q,Nt);
Field_b = zeros(q,Nt);
Spectrum = zeros(q,g);

ES = zeros(q,Nt);

%% Operators

D = -1i*2*pi*F*nf/c;
ift = exp(1i*2*pi*f*dt);

k = D'*ones(1,N);
idft = ift'*ones(1,N);
%% Initial values

Ein = ones(g,1);
Ef(:,1) = Ein;

%% loop

for r = 1:q
    l = 0;
    tic
    for n = 1:Nt
        tempf(:,1) = Ein + r1*Eb(:,end).*exp(D'*L);
        tempf(:,2:N) = (1+v*dt.*bs(:,2:N)).*Ef(:,2:N)-R*(Ef(:,2:N)-Ef(:,1:N-1));
        tempb(:,1) = r2*Ef(:,N).*exp(D'*L);
        tempb(:,2:N) = (1+v*dt.*fliplr(bs(:,2:N))).*Eb(:,2:N)-R*(Eb(:,2:N)-Eb(:,1:N-1));
        
        Ef = tempf;
        Eb = tempb;
        Es = Ef + fliplr(Eb);
        
        Eout_f = Ef(:,N)*tout;
        Eout_b = Eb(:,N)*tin;
        
        Field_f(r,n) = ift*Eout_f;
        Field_b(r,n) = ift*Eout_b;
        Et = sum((idft).*Es);
        
%         if(n*dt>l*1e-8)
%             l = l + 1;
%             subplot(311),imagesc(abs(Es))%plot(z,(abs(Et).^2))
%             subplot(312),plot(t,abs(Field(r,:)).^2)
%             subplot(313),plot(lam_s,abs(Eout))
%             pause(1e-3)
%         end
    end
    ES(r,:) = Et;
    Spectrum(r,:) = Eout';
end
figure()
t0 = 0;
for r = 1:q
plot(t0 + t,(abs(Field_f(r,n))).^2,'b',t0+t,(abs(Field_b(r,n))).^2,'g','linewidth',3)
hold on
t0 = r*T;
end

