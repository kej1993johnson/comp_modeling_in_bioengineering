% Sensitivity Analyses on drug-free data
% This script performs a sensitivity analyses on the drug-free model using
% parameter estimates from Maclean et al
%% E and M subpopulation transition model
% proliferation rates should be in range from 0.302 to 1.38 cell divisions
% per day (Maclean cites Stinson et al)
ge = 0.87; % growth rate of epithelial cells
gm = 0.32; % growth rate of mesenchymal cells
% in Maclean paper, this was varied from 0 to 0.01 { 0, 0.0001, 0.001,
% 0.0025, 0.01}
kem = 0.0025; % baseline transition rate E to M 
kme = 0.32; % transition rate from M to E


Ce_init(1) = 5e3;
Ce_init(2) = 1e3;

params = vertcat(ge, gm, kem, kme);
param_nms = {'g_{E}', 'g_{M}', 'k_{EM}','k_{ME}'} ;
tsamp = 0:0.25:6;

f = @(t,Cc) [ge*Cc(1)- kem.*Cc(1) + kme*Cc(2);  % dE/dt
             gm*Cc(2)+ kem*Cc(1) - kme*Cc(2)]; % dM/dt
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
%% Perform Sensitivity Analysis
% first do sensitivity analysis on initial conditions
figure;
for j = 1:length(params_c)
    p=params_c;
    C_i = Ce_init;
  % first find for current C_i
[ Cei] = model_EM_nodrug( C_i, params_c, tsamp);
Ni = Cei(:,1) + Cei(:,2);
% next find for increased param
    p(j)=params_c(j)*1.05;
[ Cef] = model_EM_nodrug( C_i, p, tsamp);
Nf = Cef(:,1) + Cef(:,2);
sens(:,j) = (Nf-Ni)./(params_c(j)-p(j));


plot(tsamp, sens(:,j), 'LineWidth', 2)
text(tsamp(end-2), sens(end-2,j), param_nms{j})
hold on
xlabel('time')
ylabel('Sensitivity (dN/dp)')
title('Sensitivity vectors for each parameter')

end

figure;
for i = 1:length(params_c)
    subplot(1,length(params_c),i)
    plot(tsamp, sens(:,i), 'LineWidth', 2)
    xlabel('time')
    ylabel('Sensitivity (dN/dp)')
    title([param_nms{i}])
end
    

