function[err] = fit_model_EMTpath( fit_params,outLD50, stiff,C_init, params, tsamp, beta )
%This function runs the ODEs for the mechanotransduction signaling pathway,
% given the initial concentration of each molecule and the rates

% It outputs the final concentration of each molecule at equilibrium, as
% well as the percent mesenchymal (psim) and LD50 at equlibrium
alpha = fit_params(1);
LD50base = fit_params(2);

for i =1:length(stiff)
sigma = stiff(i)
C_init(1) = sigma*alpha;
ratio = params(1);
kp = params(2);
kon = params(3);
koff = params(4);
kdep = params(5);
knuc = params(6);
kcyt = params(7);
m = params(8);


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

LD50_stiff(i)= LD50(end);
end

err = LD50_stiff-outLD50;


end

