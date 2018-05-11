function [ C, LD50, psim ] = model_EMTpath( C_init, params, tsamp, LD50base, beta )
%This function runs the ODEs for the mechanotransduction signaling pathway,
% given the initial concentration of each molecule and the rates

% It outputs the final concentration of each molecule at equilibrium, as
% well as the percent mesenchymal (psim) and LD50 at equlibrium



sigma = params(1);
alpha = params(2);
C_init(1) = sigma*alpha;
ratio = params(3);
kp = params(4);
kon = params(5);
koff = params(6);
kdep = params(7);
knuc = params(8);
kcyt = params(9);
m = params(10);


tspan =[ 0, tsamp(end)];
f = @(t,C) [0; % dKin/dt
            kon*C(3)*C(4)-koff*C(2);   % dTwG3bnd/dt
            koff*C(2)- kp*((C(3)*C(1))./(C(3)+m)) - kon*C(3)*C(4) + kdep*C(5); %dTwf/dt
            koff*C(2) - kon*C(3)*C(4); % dG3/dt
            kp*((C(3)*C(1))./(C(3)+m)) - knuc*C(5) + (kcyt*C(6)*ratio)-(kdep*C(5)); %dcTwP/dt
            knuc*C(5)*(1/ratio) - kcyt*C(6)]; %dnTwPdt

options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:6);
[t,C]=ode45(f, tsamp,C_init, options);
Kin= C(:,1);
TwG3bnd= C(:,2);
Twf = C(:,3);
G3= C(:,4);
cTwP = C(:,5);
nTwP=C(:,6);

Vnuc = 1;
Vcyt = 1/ratio;

beta = 1;
psim = ((nTwP.*Vnuc)./((Twf+TwG3bnd+cTwP).*Vcyt));
LD50 = (1 + beta*psim).*LD50base;

LD50= LD50;
psim=psim;



end

