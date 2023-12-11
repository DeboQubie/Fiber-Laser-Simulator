function [mfr,gamma]=mfd(lam,a,na)


V = 2*pi*a*na./lam;
mfr = a.*(0.65+(1.619./(V.^1.5))+(2.879./(V.^6)));
gamma = 1 - exp(-2*(a^2)./(mfr.^2));