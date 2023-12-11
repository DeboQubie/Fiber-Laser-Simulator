function dPdz = amp_2level(z,P,fiber)


%% 
c = 3e8;  % speed of light
h = 6.6256e-34; % Planck's Constant

r13 = fiber.sap*fiber.Lp/(h*c*(pi*fiber.rs^2));
r31 = fiber.sep*fiber.Lp/(h*c*(pi*fiber.rs^2));
r12 = fiber.sas.*fiber.Ls./(h*c*(pi*fiber.rs.^2));
r21 = fiber.ses.*fiber.Ls./(h*c*(pi*fiber.rs.^2));

A21 = 1/fiber.tau_f;


N1 = fiber.N0*(r31*P(1) + r21'*P(2:end) + A21)./(r13*P(1) + r31*P(1) + r12'*P(2:end) + r21'*P(2:end) + A21);
N2 = fiber.N0 - N1;

%% Gain coefficients

bp = (fiber.sep*N2 - fiber.sap*N1);
bs = (fiber.ses*N2 - fiber.sas*N1);
bn = fiber.ses*N2;

%% Differential equations

g = length(fiber.Ls);
dPdz = zeros(g+1,1);

dPdz(1,1) = bp.*P(1);                           % Pump power
dPdz(2:end,1) = bs.*P(2:end,1) + bn.*fiber.PoA;     % Signal power + ASE
% for k = 1:g
%     dPdz(k+1) = bs(k)*P(k) + bn(k)*fiber.PoA(k);
% end
