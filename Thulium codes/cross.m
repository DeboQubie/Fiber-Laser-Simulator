function out = cross(Lp,Ls,N0,r,NA)

%importing spectroscopic data from data files

h=6.626e-34;     %planck's constant in Js
c=3e8;           %speed of light in meter per sec

a=1/334.7e-6;        %1/a=lifetime of level 2 in sec
b=1/14.2e-6;       %1/b=lifetime of level 3 in sec

D=2.86e6;            %Density of Silica glass in gm/m^3
AN=6.02e23;        %Avogadro no
MM=168.9342;        %molecular mass of Tm in gm

%%
sec=importdata('Emission_crosssection_TDF.xls');    %emission cross-section at signal wavelength in square meter   % importdata("file_name");
sac=importdata('Absorption_crosssection_TDF.xls');      %absorption cross-section at signal wavelength in square meter

pindex= sac.Sheet1(:,1)==Lp;
sap=sac.Sheet1(pindex,2);

%sap=interp1(sac.Sheet1(:,1),sac.Sheet1(:,2),Lp);            %absorption cross-section at pump wavelength in square meter
sas=interp1(sac.Sheet1(:,1),sac.Sheet1(:,2),Ls);
ses=interp1(sec.Sheet1(:,1),sec.Sheet1(:,2),Ls);
%%
%calculation for signal parameters
Vs=(2*pi*r*NA)./Ls;     %V number
ws=r*(0.65+(1.619./(Vs.^1.5))+(2.879./(Vs.^6)));       %Mode field radius of signal in meters
gammas=1-exp(-2*(r^2)./(ws.^2));

%calculation for pump parameters
Vp=(2*pi*r*NA)/Lp;     %V number
wp=r*(0.65+(1.619/(Vp^1.5))+(2.879/(Vp^6)));       %Mode field radius of pump in meters
%gammap=1-exp(-2*(r^2)/(wp^2));                      %overlap factor for pump wabelength
gammap=5.917e-3;
SAP=sap*gammap;
SAS=sas.*gammas;
SES=ses.*gammas;

%calculation of cross-realaxation coefficients
wt=(N0*MM)/(AN*D);            %Dopant concentration in wt%
%K1=1.8e-22*wt/0.0575;      %K3212 cross relaxation coefficient
%K2=0;
K1=1.25e-22;
%K1=0;
K2=0.084*K1;              %K2321 cross relaxation coefficient

dL=Ls(2)-Ls(1);          %ASE channel spacing in m

PoA=2*h*(c^2)*dL./Ls.^3;            %Local noise power

out(:,1) = [Lp Ls'];
out(:,2) = [wp ws']';
out(1:306,3) = [SAP SAS']';
out(:,4) = [a SES']';
out(:,5) = [b PoA']';
out(1,6) = K1;
out(2,6) = K2;
clear sas ses sap;

