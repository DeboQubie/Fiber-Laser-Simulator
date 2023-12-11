function [N1,N2,bp,bs] = inv3(Pp,Ps,fiber)

% things it needs in fiber --> N0,K1,K2,tau21,tau32,Lp,rs,Ls,sap,sas,ses

%% 
c = 3e8;  % speed of light
h = 6.626e-34; % Planck's Constant

%% Rates

N0 = fiber.N0;
K1 = fiber.K1;
K2 = fiber.K2;

a = 1/fiber.tau21;
b = 1/fiber.tau32;

R = (Pp*fiber.Lp*fiber.sap)/(pi*(fiber.rs^2)*h*c);         % pump absorption rate, P(1)=pump power

R = R';

%fiber.sas

w21 = (fiber.Ls.*fiber.ses)./(pi.*(fiber.rs.^2)*h*c);         %Signal emission rate 
w12 = (fiber.Ls.*fiber.sas)./(pi.*(fiber.rs.^2)*h*c);         %signal absorption rate

W12 = (Ps'*w12);
W21 = (Ps'*w21);

aa = (K1*(W21+a-b).^2+K1*(W21+a-b).*(2*R+W12+b)+K2*(2*R+W12+b).^2);
bb = (2*K1*b*N0*(W21+a-b)+b*K1*N0*(2*R+W12+b)-(R+W12+N0*K1).*(W21+a-b).*(2*R+W12+b)+(W21+a).*((2*R+W12+b).^2));
cc = ((b*N0)^2*K1-b*N0*(R+W12+N0*K1).*(2*R+W12+b));

N2 = (-bb+(bb.^2-4*aa.*cc).^0.5)./(2*aa);
N1 = (N0*b+(W21+a-b).*N2)./(2*R+W12+b);

g = length(fiber.Ls);
%% Gain coefficients

%disp(["size of Ps' = ",size(Ps')]);
%disp(["size of aa = ",size(aa)]);
%disp(["size of bb = ",size(bb)]);
%disp(["size of cc = ",size(cc)]);
%disp(["size of w12 = ",size(w12)]);
%disp(["size of W12 = ",size(W12)]);
%disp(["size of W21 = ",size(W21)]);
%disp(["size of N2 = ",size(N2)]);
%disp(["size of R = ",size(R)]);
%disp(["size of K1*(W21+a-b).^2 = ",size(K1*(W21+a-b).^2)]);
%disp(["size of fiber.ses = ",size(fiber.ses)]);
%disp(["size of fiber.sas = ",size(fiber.sas)]);
%disp(["size of N1 = ",size(N1)]);
%disp(["size of K1*(W21+a-b).*(2*R+W12+b) = ",K1*(W21+a-b).*(2*R+W12+b)]);
%disp(["size of K2*(2*R+W12+b).^2 = ",K2*(2*R+W12+b).^2]);

%fiber.ses = fiber.ses';
%fiber.sas = fiber.sas';

N2 = N2';
N1 = N1';

bp = -fiber.sap*N1;
bs = fiber.ses*N2 - fiber.sas*N1;