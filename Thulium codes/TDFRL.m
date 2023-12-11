%close all
%clear all
%clc
% main function : calling function

function [Z,Ls,Spectrum,Pump_power_vs_z,Signal_power_vs_z] = TDFRL(L,N0,r,Na,alpha,Pp,lambda_p,tau21,tau32,inserting_data)

%initial conditions
% D=2.86e6;            %Density of Silica glass in gm/m^3
% AN=6.02e23;        %Avogadro no
% MM=168.9342;        %molecular mass of Tm in gm
% wt=0.007;
% N0=(wt*AN*D)/MM;

N0=N0; %2.74e26;  %Total dopant concentration in per meter cube
r=r; %5e-6;            %radius of core in meters

alpha=alpha; %5;
att=exp(-alpha/4.343);

dl=1e-9;        % wavelength spacing in m
Lp=lambda_p*1e-9; %790e-9;        %pump wavelength in m
Ls=(1850e-9:dl:2155e-9)';                    % array of signal wavelengths in m                                      
g=length(Ls);                                %no of ASE channels = 586 with spacing 1nm
% Lw=1966e-9;                                 %signal wavelength in nm
% k=find(Ls==Lw);

NA=Na; %0.15;            %Numerical Aperture of the fiber

Input = cross(Lp,Ls,N0,r,NA);

fiber_class.N0 = N0;
fiber_class.K1 = Input(1,6);
fiber_class.K2 = Input(2,6);
fiber_class.tau21 = tau21;
fiber_class.tau32 = tau32;
fiber_class.Lp = Input(1,1);
fiber_class.Ls = Input(2:end,1);
fiber_class.rs = r;
fiber_class.sap = Input(1,3);
fiber_class.sas = Input(2:end,3);
fiber_class.ses = Input(2:end,4);

fiber = Input;
fiber (3,6) = N0;
fiber (4,6) = r;

lam_p = Input(1,1);
lam_s = Input(2:end,1)';

Pp0=Pp; %4; % Initial pump power in watts
Ps0=zeros(g+1,1); % Initial signal power in watts


%initial conditions array
y0=Ps0;




%lengthspan
L=L; %2;        %length of thulium doped fiber in meters
step=100;
zspan=linspace(0,L,step);
%alpha_b=1.59;
% att_b=exp(-alpha_b/4.343);
% alpha_f=0.5;
% att_f=exp(-alpha_f/4.343);


  y0(1)=Pp0 ;
  y0(2:end)=0;
%calling ODE45
  [Z,Y]=ode45(@amp_fast,zspan,y0',fiber,fiber);

  err=8e-3;		%error value
  R=0.5;          %reflectivity of output coupler (part of signal going back into loop)
  M=1;
  i=1;

%extracting all ASE powers for each wavelength AT Z=L    
  A=Y(step,2:end);

  Pump_power_vs_z = Y(:,1);
  Signal_power_vs_z = sum(Y(:,2:end)');

  U=(1-R)*A';
plot(Ls,10*log10(U*1e3));

% drawnow
% 
  while(M>err)
      i
      M
        
%initial conditions array for backward propagation direction
   I=att*R*A';
   i0=[Pp0;I];

%calling ODE45 for forward propagation direction
   [Z,O]=ode45(@amp_fast,zspan,i0,fiber,fiber);

%extracting all ASE powers for each wavelength AT Z=L
   B=O(step,2:end);
   V=(1-R)*B';
   plot(Ls,10*log10(V*1e3));
   drawnow


%calculating error
  E=abs((V-U)./V);
  M=max(E);
  A=B;
  U=V;
  i=i+1;
  end
  

[N1,N2,bp,bs] = inv3_new(Pump_power_vs_z,Y(:,2:end)',fiber_class);

Pout=U;
Length=L;
Pump=Pp0;
Ps=max(10*log10(U*1e3));

index=find((10*log10(U*1e3))==Ps);
Wavelength=(Ls(index));
% xlswrite('Wavelength.xlsx',Ls)
% xlswrite('Spectrum_alpha_5.xlsx',10*log10(U*1e3))

plot(Ls,10*log10(U*1e3))  % Spectrum 

Spectrum = 10*log10(U*1e3);

gain = trapz(Z,bs')*4.343;

%figure 1
%plot(Ls,10*log10(U*1e3));
%plot(Pp0,Pout)
% xlswrite('Pump input power_2m_without_CR.xlsx',Pp0)
% xlswrite('Output power_2m_without_CR.xlsx',Pout)
% xlabel('Wavelength (m)','FontSize',15);
% ylabel('Laser Power (dBm)','FontSize',15);
% xlswrite('Wavlength',Wavelength)
% xlswrite('Laser power with alsph=0',10*log10(U*1e3))


if inserting_data == "y"
    
    disp(Pp);
    dlmwrite('Pp.csv',Pump_power_vs_z);
    dlmwrite('Ps.csv',Y(:,2:end));
    dlmwrite('N1.csv',N1);
    dlmwrite('N2.csv',N2);
    dlmwrite('bp.csv',bp);
    dlmwrite('bs.csv',bs);
    dlmwrite('lams.csv',Ls);

    disp("Initial conditions written into the .csv files - Closed");

else 

    figure()
    subplot(221),plot(Z,Pump_power_vs_z),xlabel('Length (m)'),ylabel('Power (mW)')%,Z,Power*1e3),title('Power vs Length'),legend('Pump','Signal'),grid on
    subplot(222),plot(Z,Signal_power_vs_z),xlabel('Length (m)'),ylabel('Power (mW)')
    subplot(223),plot(Ls*1e9,gain),xlabel('Wavelength (nm)'),ylabel('Gain (dB)')
    subplot(224),plot(Ls,10*log10(U*1e3)),xlabel('Wavelength (m)'),ylabel('Power (dBm)')%,lam_s,Noise)
    pause(5e-3)

end 

%drawnow




