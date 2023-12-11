function dPdz = amp3(z,P,fiber)


%% 
c = 3e8;  % speed of light
h = 6.626e-34; % Planck's Constant

%% Rates

N0 = fiber.N0;
K1 = fiber.K1;
K2 = fiber.K2;

a = 1/fiber.tau21;
b = 1/fiber.tau32;

R = (P(1)*fiber.Lp*fiber.sap)/(pi*(fiber.rs^2)*h*c);         % pump absorption rate, P(1)=pump power

w21 = (fiber.Ls.*fiber.ses)./(pi.*(fiber.rs.^2)*h*c);         %Signal emission rate 
w12 = (fiber.Ls.*fiber.sas)./(pi.*(fiber.rs.^2)*h*c);         %signal absorption rate

W12 = P(2:end)'*w12;
W21 = P(2:end)'*w21;

aa = (K1*(W21+a-b).^2+K1*(W21+a-b).*(2*R+W12+b)+K2*(2*R+W12+b).^2);
bb = (2*K1*b*N0*(W21+a-b)+b*K1*N0*(2*R+W12+b)-(R+W12+N0*K1).*(W21+a-b).*(2*R+W12+b)+(W21+a).*((2*R+W12+b).^2));
cc = ((b*N0)^2*K1-b*N0*(R+W12+N0*K1).*(2*R+W12+b));

N2 = (-bb+(bb.^2-4*aa.*cc).^0.5)./(2*aa);
N1 = (N0*b+(W21+a-b).*N2)./(2*R+W12+b);

g = length(fiber.Ls);
%% Gain coefficients

bp = -fiber.sap*N1;
bs = fiber.ses*N2 - fiber.sas*N1;
bn = fiber.ses*N2;

%% Differential equations

dPdz = zeros(g+1,1);

dPdz(1,1) = bp.*P(1);                           % Pump power
dPdz(2:end,1) = bs.*P(2:end,1) + bn.*fiber.PoA;     % Signal power + ASE

% 
% initial=y(2:end);
% 
% output = (bs.*initial)+ (bn.*PoA);   %signal+ASE power
% dydt(2:end,1)=output;