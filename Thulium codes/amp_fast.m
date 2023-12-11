%called function
function dydt=amp_fast(z,y,fiber)
  
c=3e8;  % speed of light
h=6.626e-34; % Planck's Constant
%%
r=fiber(4,6);         %core radius
Lp=fiber(1,1);      %pump wavelength

wp=fiber(1,2);      %mode field radius at pump wavelength

N0=fiber(3,6);        %total dopant concentration
sap=fiber(1,3);     %abs cross section at pump wavelength

Ls=fiber(2:end,1);      %Laser wavelength array

ws=fiber(2:end,2);          %mode field radius at signal wavelengths

sas=fiber(2:end,3);         %absorption cross-section at signal wavelengths
ses=fiber(2:end,4);         %emission cross-section at signal wavelengths

PoA=fiber(2:end,5);             %Local noise power

a=fiber(1,4);        %1/a=lifetime of level 2 in sec
b=fiber(1,5);       %1/b=lifetime of level 3 in sec

K1=fiber(1,6);  %cross-realaxation coefficient 1
K2=fiber(2,6);  %cross-realaxation coefficient 2

g=length(Ls);
%%
%pump and signal absorption rates
R=(y(1)*Lp*sap)/(pi*(r^2)*h*c);         %pump absorption rate, y(1)=pump power
                                         
w21=(Ls.*ses)./(pi.*(r.^2)*h*c);         %Signal emission rate 
w12=(Ls.*sas)./(pi.*(r.^2)*h*c);         %signal absorption rate

%adding all the signal absorption and emmision for all wavelengths to
W12 = y(2:end)'*w12;
W21 = y(2:end)'*w21;

%Dopant concentrations
% N2 = N0*b*(R+W12)/(R*b+R*a+a*b+R*W21+b*W21+W12*b);
% N1 = b*(N0-N2)/(R+b);

aa=(K1*(W21+a-b).^2+K1*(W21+a-b).*(2*R+W12+b)+K2*(2*R+W12+b).^2);
bb=(2*K1*b*N0*(W21+a-b)+b*K1*N0*(2*R+W12+b)-(R+W12+N0*K1).*(W21+a-b).*(2*R+W12+b)+(W21+a).*((2*R+W12+b).^2));
cc=((b*N0)^2*K1-b*N0*(R+W12+N0*K1).*(2*R+W12+b));

N2=(-bb+(bb.^2-4*aa.*cc).^0.5)./(2*aa);
N1=(N0*b+(W21+a-b).*N2)./(2*R+W12+b);

%Differential equations
bp=-sap*N1;
bs=ses*N2-sas*N1;
bn = ses*N2;
dydt = zeros(g+1,1);
dydt(1,1)=bp.*y(1);               %pump power

initial=y(2:end);
output = (bs.*initial)+ (bn.*PoA);   %signal+ASE power
dydt(2:end,1)=output;






