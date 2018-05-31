 %This code performs simulations of the YAP/TAZ mechanotransduction
% signaling pathway and its effect on the subpopulation composition of E
% and M cells.

% We first hypothesize that YAP/TAZ nuclear localization directly modulates
%  downstream markers of E and M cell, and the ratio of these markers 
% determines the  rate of E to M transition from the Maclean et al 2014 paper.
% YAP/TAZ nuclear localization changes in response to a mechnotransduction
% signaling pathway.
close all; clear all; clc

% Subcellular model & population level model combines
% 1. YAPnuc
% 2. YAPcyt
% 3. Vim
% 4. Ecad
% 5. E cells
% 6. M cells

%% Start with subcellular model
tsamp = 0:1:144;
% input mechanical function (for static this should be just a step function
% or constant, but eventually will make stiffness go up and down)
sigma = 2000;%*ones(length(tsamp), 1); % mechanical stress function 
% rate parameters for subscellular model
kcn = 0.0001; % rate constant determining effect of stiffness on rate of transport of YAP to nucleus
knc = 0.8; % rate of transport of YAP from the nucleus to the cytoplasm
Vnc = 0.2; % ratio of volume of nucleus to the cytoplasm
kv= 0.002; % rate of additional Vim production due to nuclear localization of YAP/TAZ
dv = 0.0001; % rate of degradation of Vim
pv = 0.0001; % rate of production of Vim
kE= 0.1; % rate of additional Ecad degradation due to nuclear localization of YAP/TAZ
dE = 0.0001; % rate of degradation of Ecad
pE = 0.0001; % rate of production of Ecad
% rate parameters for population level model
ge = 0.037; % growth rate of epithelial cells
kem = 0.4; % baseline transition rate E to M
kme = 0.32; % transition rate from M to E
gm = 0.032; % growth rate of mesenchymal cells
carcap = 1e4;

C_init(1) = 0.1; % initial fraction(concentration?) YAP in nucleus
C_init(2) = 0.9; % initial fraction (concentration?) YAP in cytoplasm
C_init(3)= 0.01; % inital concentration of Vimentin (very low)
C_init(4) = 0.8; % inital concentration of E-cadherin (high)
C_init(5) = 800; % E cells
C_init(6) = 200; % M cells

params_ = vertcat(sigma, kcn, knc, Vnc, kv, dv,pv, kE,dE, pE, ge, kem, kme, gm, carcap);
param_nms = {'\sigma', 'k_{cyt-nuc}', 'k_{nuc-cyt}','V_{nuc-cyt}','k_{p}', 'k_{Vim}', 'deg_{Vim}','prod_{Vim}', 'k_{Ecad}', 'def_{Ecad}', 'prod_{Ecad}', 'g_{E}', 'k_{EM}', 'k_{ME}', 'g_{M}'} ;
tsamp = 0:1:144;

f = @(t,C) [kcn*sigma.*C(2).*((Vnc)^-1) - knc.*C(1).*Vnc; % dYAPnuc/dt
            -kcn*sigma.*C(2).*((Vnc)^-1) + knc.*C(1).*Vnc;   % dYAPcyt/dt
            kv.*(C(1)) - dv + pv; % dVim/dt
            -kE.*(C(1)) - dE + pE; % dEcad/dt
            ge*(1-((C(5)+C(6))./carcap))*C(5)- kem.*(C(1)./(C(1)+C(2))).*C(5) + kme*C(6);  % dE/dt
             gm*(1-((C(5)+C(6))./carcap))*C(6)+ kem*(C(1)./(C(1)+C(2))).*C(5) - kme*C(6)]; % dM/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:6);
[t,C]=ode45(f, tsamp,C_init, options);
YAPnuc= C(:,1);
YAPcyt= C(:,2);
Vim = C(:,3);
Ecad = C(:,4);
E = C(:,5);
M = C(:,6);
names = {'[YAP_{nuc}]', '[YAP_{cyt}]','[Vim]', '[Ecad]', 'E cells', 'M cells'};

figure;
for i = 1:length(C_init)
subplot(2,3,i)
plot(tsamp, C(:,i), 'LineWidth', 2)
xlabel('time')
ylabel('concentration')
title(names(i))
xlim([0, tsamp(end)])
end

figure
plot(tsamp, YAPnuc./(YAPcyt+YAPnuc), 'LineWidth',2)
hold on
plot(tsamp, Vim, 'LineWidth',2)
plot(tsamp, Ecad, 'LineWidth',2)
xlabel('time')
ylabel('fraction YAP in nucleus')
legend('fraction YAP_{nuc}', '[Vim]', '[Ecad]')
title('Time course of EMT drivers and markers')
legend boxoff

figure
plot(tsamp, E, 'b-','LineWidth',2)
hold on
plot(tsamp, M, 'r-','LineWidth',2)
xlabel ('time')
ylabel('number of cells')
title('E and M cells over time as YAP nuclear localization increases soft conditions')
legend ('E cells', 'M cells')
legend boxoff

fres = M(end)/(E(end)+M(end))

