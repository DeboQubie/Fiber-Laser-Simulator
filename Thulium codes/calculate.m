

function [N1,N2,N3,CR] = calculate(Pp,Ps)

%Parameters
c=3e8;  % speed of light
h=6.625e-34; % Planck's Constant
%%
Input = dlmread('input.csv');
fiber = dlmread('fiber.csv');

r=fiber(2);         %core radius
Lp=Input(1,1);      %pump wavelength

wp=Input(1,2);      %mode field radius at pump wavelength

N0=fiber(1);        %total dopant concentration
sap=Input(1,3);     %abs cross section at pump wavelength

Ls=Input(2:end,1);

ws=Input(2:end,2);

sas=Input(2:end,3);
ses=Input(2:end,4);

PoA=Input(2:end,5);

a=Input(1,4);
b=Input(1,5);

K1=Input(1,6);
K2=Input(2,6);

g=length(Ls);
clear Input fiber;


%pump and signal absorption rates
R=(Pp*Lp*sap)/(pi*(r)^2*h*c);         %pump absorption rate
                                                                                
w21=(Ls.*ses)./(pi.*(r.^2)*h*c);         %Signal emission rate
                                                                                   
w12=(Ls.*sas)./(pi.*(r.^2)*h*c);         %signal absorption rate

%adding all the signal absorption and emmision for all wavelengths
g21=Ps*w21;
t12=Ps*w12;

W12 = t12';
W21 = g21';

%Dopant concentrations
% N2 = N0*b*(R+W12)/(R*b+R*a+a*b+R*W21+b*W21+W12*b);
% N1 = b*(N0-N2)/(R+b);
%N2=((-(2*N0*b*K1*(W21+a-b)+K1*N0*b*(W12+2*R+b)+(-R-W12-N0*K1).*(W21+a-b).*(2*R+W12+b)+(W21+a).*(2*R+W12+b).^2))+((2*N0*b*K1*(W21+a-b)+K1*N0*b*(W12+2*R+b)+(-R-W12-N0*K1).*(W21+a-b).*(2*R+W12+b)+(W21+a).*(2*R+W12+b).^2).^2-4*(K1*(W21+a-b).^2+K1*(W21+a-b).*(2*R+W12+b)+K2*(2*R+W12+b).^2).*(K1*(N0*b)^2+(-R-W12-N0*K1).*(2*R+W12+b)*N0*b)).^0.5)./(2*(K1*(W21+a-b).^2+K1*(W21+a-b).*(2*R+W12+b)+K2*(2*R+W12+b).^2));
%N1=(N0*b+(W21+a-b).*N2)./(2*R+W12+b);
 
aa=(K1*(W21+a-b).^2+K1*(W21+a-b).*(2*R+W12+b)+K2*(2*R+W12+b).^2);
bb=(2*K1*b*N0*(W21+a-b)+b*K1*N0*(2*R+W12+b)-(R+W12+N0*K1).*(W21+a-b).*(2*R+W12+b)+(W21+a).*((2*R+W12+b).^2));
cc=((b*N0)^2*K1-b*N0*(R+W12+N0*K1).*(2*R+W12+b));

N2=(-bb+(bb.^2-4*aa.*cc).^0.5)./(2*aa);
N1=(N0*b+(W21+a-b).*N2)./(2*R+W12+b);

N3=N0-(N1+N2);
CR=N0.*N1.*K1-(N1).^2*K1-N2.*N1.*K1-(N2).^2.*K2;
bp=-sap*N1;
bs=ses*N2-sas*N1;