close all
clear all
clc
%% channels

EDFA_L = 3.5;                                   % length of EDFA
EDFA_w = 4.4e-6/2;                                  % mode field radius
NA = 0.29;                                      % numerical aperture
EDFA_N0 = 2.0e25;

dlmwrite('fiber.csv',[EDFA_N0,EDFA_w]);
fiberp = [EDFA_N0,EDFA_w];

wdm_p = 1.55;
iso_l = 0.9;
cou_l = 0;
% wdm_l = 1.6;
% aom_l = 0;%3.35;
% c_l = 0;%7;

% Cav_Loss = iso_l + cou_l + wdm_l + aom_l + c_l;       %Cavity Loss in dB
% loss = exp(-Cav_Loss/4.343);
rout = exp(-(iso_l + cou_l)/4.343);

lam_p = 976;      % Pump wavelength in nm

lam1 = 1501;        % ASE start wavelength in nm

lam2 = 1600;        % ASE end wavelength in nm
dl = 0.1;             % ASE wavelength separation in nm

lam_s = lam1:dl:lam2;
lam = [lam_p lam_s];
% dlmwrite('lam.csv',[lam_p lam_s])

Input = spect(lam,EDFA_N0,EDFA_w,NA);
dlmwrite('input.csv',Input);
%%

lam_p = Input(1,1);
lam_s = Input(2:end,1)';

Noise = Input(2:end,5)';

% q = 1000;
err = 1;eps = 1e-3;
r=0;Pp = 0;Ps = 0;N1 = 0;
% Po = zeros(1,5000);
%% 

N = 466;
dz = EDFA_L/N;

z = (1:N)'*dz;

Pump = 1e-3*input('Pump power (mW) = ');
Signal = 00e-6;
sig_lam = 1550e-9;
P_in = Pump*exp(-wdm_p/4.343);
S_in = zeros(1,length(lam_s));
S_in(abs(lam_s-sig_lam)<(1e-9*dl/2)) = Signal;
disp(['Pump coupled in EDF = ',num2str(P_in*1e3),' mW'])
disp(['Signal wavelength = ',num2str(sig_lam*1e9),' nm'])
disp(['Input Power = ',num2str(Signal*1e3),' mW'])
% while (err>=eps)
%     r=r+1;

%     tic
%     PP = Pp;PS = Ps;n1 = N1;
    [Z,Sol] = ode45(@CW_nlam,z,[P_in S_in],Input,Input,fiberp);
    
    Z=Z';
    
    Pp = Sol(:,1)';
    Ps = Sol(:,2:end)';
    Power = sum(Ps);
    P_end = Ps(:,end)';
    Spectrum = P_end*rout;
%     S_in = loss*P_end;
    [N1,bp,bs] = inversion_n(Pp,Ps,Input,fiberp);
%     err = max([max(abs(Pp-PP)./PP),max(max(abs(Ps-PS)./PS))]);
%     disp(r),disp(err)
    Po = sum(Spectrum);
    figure()
    subplot(221),plot(Z,Pp*1e3),xlabel('Length (m)'),ylabel('Pump power (mW)')%,Z,Power*1e3),title('Power vs Length'),legend('Pump','Signal'),grid on
    subplot(222),plot(Z,[Power;Ps(abs(lam_s-sig_lam)<(1e-9*dl/2),:)]*1e3),xlabel('Length (m)'),ylabel('Signal power (mW)')
    subplot(223),plot(Z,[N1;1-N1]),axis([0 EDFA_L 0 1]),xlabel('Length (m)'),ylabel('N_1/N_0, N_2/N_0')%,title('Excited State Population'),grid on
    subplot(224),plot(lam_s,10*log10(Spectrum*1e3)),xlabel('Wavelength (m)'),ylabel('Power (mW)')%,lam_s,Noise)
%     pause(5e-3)
%     toc
% end
disp(['Output Power = ',num2str(Po*1e3),' mW'])
disp(['Residual Pump = ',num2str(Pp(end)*1e3),' mW'])
% subplot(221),plot(Z,Pp*1e3),xlabel('Length (m)'),ylabel('Power (mW)'),title('Pump Power vs Length'),grid on
% subplot(222),plot(Z,Power*1e3),xlabel('Length (m)'),ylabel('Power (mW)'),title('Signal Power vs Length'),grid on
% subplot(223),plot(Z,1-N1),axis([0 2 0 1]),xlabel('Length (m)'),ylabel('N_2/N_0'),title('Excited State Population'),grid on
% subplot(224),plot(lam_s*1e9,10*log10(Spectrum)),xlabel('Wavelength (nm)'),ylabel('Power (dBm)'),title('Spectrum'),grid on
% figure(),plot((1:r),Po(1:r)),ylabel('Power (W)'),title('Power'),grid on
gain = trapz(z,bs')*4.343 ;
figure(),plot(lam_s*1e9,gain),xlabel('Wavelength (m)'),ylabel('Gain (dB)')
disp(['Signal Gain = ',num2str(gain(abs(lam_s-sig_lam)<1e-16)),' dB'])
gs = 4.343*bs(abs(lam_s-sig_lam)<1e-15,:);
figure(),plot(Z,gs)
% dlmwrite('Pp.csv',Pp);
% dlmwrite('Ps.csv',Ps);
% dlmwrite('N1.csv',N1);
disp(['Pump absorption = ',num2str((-10*log10((Pp(end))/P_in))/EDFA_L),' dB/m'])
% dlmwrite('Pp.csv',Pp);
% dlmwrite('Ps.csv',Ps);
% dlmwrite('N1.csv',N1);
% dlmwrite('bp.csv',bp);
% dlmwrite('bs.csv',bs);
% dlmwrite('lams.csv',lam_s);
delete input.csv fiber.csv