function [N1,bp,bs] = inv_2level(Pp,Ps,fiber)


%% 
c = 3e8;  % speed of light
h = 6.626e-34; % Planck's Constant

r13 = fiber.sap*fiber.Lp/(h*c*(pi*fiber.rs^2));
r31 = fiber.sep*fiber.Lp/(h*c*(pi*fiber.rs^2));
r12 = fiber.sas.*fiber.Ls./(h*c*(pi*fiber.rs.^2));
r21 = fiber.ses.*fiber.Ls./(h*c*(pi*fiber.rs.^2));

A21 = 1/fiber.tau_f;

N1 = fiber.N0*(r31*Pp + r21'*Ps + A21)./(r13*Pp + r31*Pp + r12'*Ps + r21'*Ps + A21);
N2 = fiber.N0 - N1;

%% Gain coefficients

bp = (fiber.sep*N2 - fiber.sap*N1);
bs = (fiber.ses*N2 - fiber.sas*N1);
