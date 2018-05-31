% This code performs simulations of the YAP/TAZ mechanotransduction
% signaling pathway and its effect on the subpopulation composition of E
% and M cells.

% We first hypothesize that YAP/TAZ nuclear localization directly modulates
%  downstream markers of E and M cell, and the ratio of these markers 
% determines the  rate of E to M transition from the Maclean et al 2014 paper.
% YAP/TAZ nuclear localization changes in response to a mechnotransduction
% signaling pathway.

% Subcellular model
% 1. YAPnuc
% 2. YAPcyt
% 3. Vim
% 4. Ecad

% Population Level model
% 1. E cells
% 2. M cells

%% Start with subcellular model
sigma = 200; % mechanical stress
kcn = 0.001; % rate constant determining effect of stiffness on rate of transport of YAP to nucleus
knc = 0.8; % rate of transport of YAP from the nucleus to the cytoplasm
Vnc = 0.2; % ratio of volume of nucleus to the cytoplasm
kv= 0.002; % rate of additional Vim production due to nuclear localization of YAP/TAZ
dv = 0.0001; % rate of degradation of Vim
pv = 0.0; % rate of production of Vim
kE= 0.1; % rate of additional Ecad degradation due to nuclear localization of YAP/TAZ
dE = 0.0; % rate of degradation of Ecad
pE = 0.0001; % rate of production of Ecad

C_init(1) = 0.1; % initial fraction(concentration?) YAP in nucleus
C_init(2) = 0.9; % initial fraction (concentration?) YAP in cytoplasm
C_init(3)= 0.01; % inital concentration of Vimentin (very low)
C_init(4) = 0.8; % inital concentration of E-cadherin (high)

params = vertcat(sigma, kcn, knc, Vnc, kv, dv,pv, kE,dE, pE);
param_nms = {'\sigma', 'k_{cyt-nuc}', 'k_{nuc-cyt}','V_{nuc-cyt}','k_{p}', 'k_{Vim}', 'deg_{Vim}','prod_{Vim}', 'k_{Ecad}', 'def_{Ecad}', 'prod_{Ecad}'} ;
tsamp = 0:1:144;

f = @(t,C) [kcn*sigma.*C(2).*((Vnc)^-1) - knc.*C(1).*Vnc; % dYAPnuc/dt
            -kcn*sigma.*C(2).*((Vnc)^-1) + knc.*C(1).*Vnc;   % dYAPcyt/dt
            kv.*(C(1)) - dv + pv;
            -kE.*(C(1)) - dE + pE];
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:4);
[t,C]=ode45(f, tsamp,C_init, options);
YAPnuc= C(:,1);
YAPcyt= C(:,2);
Vim = C(:,3);
Ecad = C(:,4);
names = {'[YAP_{nuc}]', '[YAP_{cyt}]','[Vim]', '[Ecad]'};

figure;
for i = 1:length(C_init)
subplot(2,2,i)
plot(tsamp, C(:,i), 'LineWidth', 2)
xlabel('time')
ylabel('concentration')
title(names(i))
xlim([0, tsamp(end)])
end

figure
plot(tsamp, YAPnuc./YAPcyt, 'LineWidth',2)
hold on
plot(tsamp, Vim./Ecad, 'LineWidth',2)
xlabel('time')
ylabel('fraction YAP in nucleus')

%% E and M subpopulation transition model
Vf = Vim(end)./Vim(1);
ge = 0.57; % growth rate of epithelial cells
beta = 0.00025; % baseline transition rate E to M
ro = 0.32; % transition rate from M to E
gm = 0.32; % growth rate of mesenchymal cells
carcap = 1e6;

Ce_init(1) = 1;
Ce_init(2) = 1;

params_c = vertcat(ge, beta, ro, gm);
param_nms = {'g_{E}', '\beta', '\ro', 'g_{M}'} ;
tsamp = 0:1:72;

f = @(t,Cc) [ge*(1-(Cc(1)./carcap))*Cc(1)- beta.*Cc(1) + ro*Cc(2);  % dE/dt
             gm*Cc(2)+ beta*Vf.*Cc(1) - ro*Cc(2)]; % dM/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:2);
[t,Cc]=ode45(f, tsamp,Ce_init, options);
namesc={'E', 'M'};
E= Cc(:,1);
M = Cc(:,2);

figure;
plot(tsamp, Cc(:,1), 'b', 'LineWidth', 2)
hold on
plot(tsamp, Cc(:,2), 'r', 'LineWidth',2)
xlabel ('time')
ylabel('cells')
legend ('E cells', 'M cells')
legend boxoff
title('E and M cell trajectories for soft matrices')

fres = M(end)./(E(end)+M(end))


