% This code performs simulations of the mechanotransduction signaling
% pathway.

close all; clear all; clc
% molecules to be modeled:
% In ODEs
% 1. Kinase (Kin)
% 2. Cytoplasmicly bound Tw:G3 (TwG3bnd)
% 3. Free cytoplasmic Twist (Twf)
% 4. Unbound G3 (G3)
% 5. Cytoplasmic phosphorylated Twist (cTwP)
% 6. Nucleus phosphorylated Twist (nTwP)

% As a result of ODEs
% 7. Mesenchymal-ness (proportional to TwPnuc/Total Tw)
% 8. LD50 proportional to mesenchymalness of cells
%%

% Assume at onset, low cell stress (sigma), low TW:P nucleus

sigma = 200; % mechanical stress
alpha = 0.0001; % proportionality constant
C_init(1)=sigma*alpha; % Kinase (t=0)
C_init(2)=0.01; % TwG3bnd (t=0)
C_init(3)= 0.5; % Twf (t=0)
C_init(4)= 0.8; % G3(t=0)
C_init(5) = 0.01; % cTwP (t=0)
C_init(6)= 0.01; % nTwP

ratio = 0.2;% Vnuc to Vcyt

 % average mechanical stress on cell due to stiffness of substrate
             % this will scale with stiffness of substrate (soft  is 200
             % Pa, hard will be 2000)
kp = 0.5; % phosphorylation rate, assume is direction proportional to average stress
kon = 0.7; % rate Twist binds to G3P2 anchored to cytoplasm
koff = 0.1; % rate Twist ubninds from G3P2
kdep= 0.01; % rate of dephosphorylation of Twist
knuc= 0.4; % rate of diffusion of phosporylated Twist into the nucleus
kcyt = 0.4; % rate of diffusion of phosphorylated Twist out of nucleus
m=C_init(3); % michaelis-menten constant, usualy set equal to concentration of substrate
params = vertcat(sigma, alpha, ratio,kp, kon, koff, kdep, knuc, kcyt, m);
param_nms = {'\sigma', '\alpha', 'V_{ratio}','k_{p}', 'k_{on}', 'k_{off}', 'k_{dep}', 'k_{nuc}', 'k_{cyt}', 'm'};
tsamp = 0:1:144; % hours? need to look at hunter's data
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

beta = 1;
LD50base = 30;
Vnuc = 1;
Vcyt = 5;
psim = ((nTwP.*Vnuc)./((Twf+TwG3bnd+cTwP).*Vcyt));
LD50 = (beta*psim +1).*LD50base;

LD50eq= LD50(end);
psimeq=psim(end);

Cend= C(end, :);

names = { '[K]', '[Tw:G3]', '[Tw]', '[G3]', '[TwP]','[TwP_{nuc}]'};
%%
figure;
for i = 1:length(C_init)
subplot(2,3,i)
plot(tsamp, C(:,i), 'LineWidth', 2)
xlabel('time')
ylabel('concentration')
title(names(i))
xlim([0, tspan(end)])
end

figure;
subplot(1,3,1)
hold off
plot(tsamp, C(:,6), 'LineWidth',2)
xlabel('time')
ylabel('concentration')
title('Nuclear Twist-P')
xlim([0, tspan(end)])
subplot(1,3,2)
hold off
plot(tsamp, psim, 'LineWidth',2)
xlabel('time')
ylabel('Mesenchymal-ness')
title('Mesenchymal-ness')
xlim([0, tspan(end)])
subplot(1,3,3)
hold off
plot(tsamp, LD50, 'LineWidth',2)
xlabel('time')
ylabel('LD50')
title('LD50')
xlim([0, tspan(end)])
%% Test out function
% just check to make sure you get the same results
[ Ctest, LD50test, psimtest ] = model_EMTpath( C_init, params, tsamp, LD50base, beta );
%% What is relationship between stiffness and LD50?

stiffvec = linspace(200, 3000, 25)
p = params;
for j = 1:length(stiffvec)
    p(1) = stiffvec(j)
    [ C, LD50stiff, psimstiff ] = model_EMTpath( C_init, p, tsamp, LD50base, beta );
    LD50endst(j) = LD50stiff(end);
end

figure;
plot(stiffvec, LD50endst, 'LineWidth', 2,'color','r')
xlabel('matrix stiffness \sigma (Pa)')
ylabel('LD50 (\muM) at 6 days')
xlim([ 200 3000])
title ('LD50 as a function of stiffness of matrix')
%% Exp 1 Use function to test effect of G3BP2 knockdown in soft and stiff conditions
G3BP2base = C_init(4);
G3BP2vec = linspace(0,G3BP2base, 5);
sigmavec = [200, 3000];
C_i = C_init;
p = params;
p(1) = sigmavec(1);
for i= 1:length(G3BP2vec)
    C_i(4) = G3BP2vec(i);
    [ C1, LD501(:,i), psim1(:,i) ] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
    C_nTwp(:,i) = C1(:,6);
    
end
p2 = params;
p2(1) = sigmavec(2);
for i= 1:length(G3BP2vec)
    C_i(4) = G3BP2vec(i);
    [ C2, LD502(:,i), psim2(:,i) ] = model_EMTpath( C_i, p2, tsamp, LD50base, beta );
    C_nTwp2(:,i) = C2(:,6);
end
%%
figure; hold 
for i =1:length(G3BP2vec)
    
    subplot(2,2,1)
    hold on
    plot(tsamp, C_nTwp(:,i), 'LineWidth',2)
    text(tsamp(end-3), C_nTwp(end-5,i), ['[G3BP2]_{0}=', num2str((G3BP2vec(i)))],'HorizontalAlignment','right','VerticalAlignment','bottom')
    xlabel('time')
    ylabel('[TwP]_{nuc}')
    title('Soft \sigma = 200 Pa')
    xlim([0, tspan(end)])
end
    subplot(2,2,2)
    hold on
    bar(G3BP2vec, psim1(end,:), 'LineWidth',2)
    %text(tsamp(end-3), psim1(end-3,i), ['[G3BP2]_{0}=', num2str((G3BP2vec(i)))],'HorizontalAlignment','right','VerticalAlignment','bottom')
    xlabel('[G3BP2]_{0}')
    ylabel('mesenchymalness (\psi)')
    title('Soft \sigma = 200 Pa')
    xlim([ G3BP2vec(1)-.1, G3BP2vec(end)+.1])
    
for i =1:length(G3BP2vec)
    
    subplot(2,2,3)
    hold on
    plot(tsamp, C_nTwp2(:,i), 'LineWidth',2)
    %text(tsamp(end-3), C_nTwp2(end-5,i), ['[G3BP2]_{0}=', num2str((G3BP2vec(i)))],'HorizontalAlignment','right','VerticalAlignment','bottom')
    xlabel('time')
    ylabel('[TwP]_{nuc}')
    title('Hard \sigma = 3000 Pa')
    xlim([0, tspan(end)])
end    

    subplot(2,2,4)
    hold on
    bar(G3BP2vec, psim2(end,:), 'LineWidth',2)
    %text(tsamp(end-3), psim1(end-3,i), ['[G3BP2]_{0}=', num2str((G3BP2vec(i)))],'HorizontalAlignment','right','VerticalAlignment','bottom')
    xlabel('[G3BP2]_{0}')
    ylabel('mesenchymalness (\psi)')
    title('Hard \sigma = 3000 Pa')
    xlim([ G3BP2vec(1)-.1, G3BP2vec(end)+.1])
%% Next perform similar analysis to model effect on TWIST1 levels on mesenchymalness

Twfbase = C_init(3);
Twfvec = linspace(0,Twfbase, 5);
sigmavec = [750, 1000];
C_i = C_init;
p = params;
p(1) = sigmavec(1);
for i= 1:length(Twfvec)
    C_i(3) = Twfvec(i);
    [ C3, LD503(:,i), psim3(:,i) ] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
    C_nTwp3(:,i) = C3(:,6);
    
end
p2 = params;
p2(1) = sigmavec(2);
for i= 1:length(Twfvec)
    C_i(4) = Twfvec(i);
    [ C4, LD504(:,i), psim4(:,i) ] = model_EMTpath( C_i, p2, tsamp, LD50base, beta );
    C_nTwp4(:,i) = C4(:,6);
end
%%
figure;  
for i =1:length(Twfvec)
    
    subplot(1,2,1)
    hold on
    plot(tsamp, C_nTwp3(:,i), 'LineWidth',2)
    text(tsamp(end-3), C_nTwp3(end-5,i), ['[TWIST1]_{0}=', num2str((Twfvec(i)))],'HorizontalAlignment','right','VerticalAlignment','bottom')
    xlabel('time')
    ylabel('[TwP]_{nuc}')
    title('Time course of TWIST1 nuclear localization')
    xlim([0, tspan(end)])
end
    subplot(1,2,2)
    hold on
    bar(Twfvec, C_nTwp3(end,:), 'LineWidth',2)
    %text(tsamp(end-3), psim1(end-3,i), ['[G3BP2]_{0}=', num2str((G3BP2vec(i)))],'HorizontalAlignment','right','VerticalAlignment','bottom')
    xlabel('[TWIST1]_{0}')
    ylabel('[TwP]_{nuc}')
    title('TWIST1 in nucleus at end')
    xlim([ Twfvec(1)-.1, Twfvec(end)+.1])
%% Sensitivity analysis
% first do sensitivity analysis on initial conditions
p=params;
for j = 2:length(C_init)
    C_i = C_init;
  % first find for current C_i
[ Cs, LD50curr, psims] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
% next find for incerased param
    C_i(j)=C_init(j)*1.01;
    deltaC=0.01*C_init(j);
[ Cs, LD50sens, psims] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
sens_C(j-1) = norm(LD50sens-LD50curr)/deltaC;
end

figure;
subplot(1,2,1)
barh(1:1:length(C_i)-1, ((sens_C)./max(sens_C)))
set(gca,'yticklabel',names(2:end))
ylabel('Protein')
xlabel('Sensitivity Score')
title('A  Sensitivity Analysis on LD50 in Initial Concentrations')
C_i = C_init;
for j = 2:length(params)
    p = params;
  % first find for current p
[ Cs, LD50curr, psims] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
% next find for incerased param
    p(j)=params(j)*1.01;
    deltap=0.01*params(j);
[ Cs, LD50sens, psims] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
sens_p(j-1) = norm(LD50sens-LD50curr)/deltap;
end
subplot(1,2,2)
barh(1:1:length(p)-1, ((sens_p)./max(sens_p)))
set(gca,'yticklabel',param_nms(2:end))
ylabel('Parameter')
xlabel('Sensitivity Score')
title('B   Sensitivity Analysis on LD50 in Parameters')


%% Now vary a few parameters and see how it affects
alphavec = linspace(0, 10*alpha, 100);
kdepvec = linspace(0, 10*kdep, 100);

C_i = C_init;
p = params;

for j = 1:length(alphavec)
    p(2) = alphavec(j);
for i = 1:length(kdepvec)
    p(7) = kdepvec(i);
    [ Cendtest, LD50eq, psimeqtest] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
LD50eqtestmat(j,i) = LD50eq(end);
end 
end
%%
imagesc(LD50eqtestmat)
colorbar
ylabel ('\alpha')
xlabel('k_{dep}')
% set(gca,'yticklabel',alphavec)
% set(gca,'xticklabel',kdepvec)
title('Map of parameter space \alpha versus kdep on LD50')
%% Model calibration to Hunter's data
stiff(:,1) = [300, 2000, 3000];
outLD50(:,1) = [6, 12, 16.6];
%outLD50(:,1) = [9.35, 12.09, 20.89];
const_params = params(3:end);


LB = [0 0];  % Lower Bounds
UB = [Inf Inf]; % Upper Bounds
params0 = [params(2) 8];% Initial Guess for alpha and LD50 base
options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
% test error function
%err = fit_model_EMT( params0,outLD50, stiff,C_init, const_params, tsamp, beta )
% Now fit to only two of the stiffnesses and LD50s
[fit_params, resnorm, reslsq]= lsqnonlin(@fit_model_EMT, params0, LB, UB, options, outLD50(1:2), stiff(1:2),C_init, const_params, tsamp, beta );
Rsq = 1- (sum((reslsq.^2))./(sum((mean(outLD50(1:2))-outLD50(1:2)).^2)))

% Model calibration
p = params;
p(2) = fit_params(1); % replace alpha with model fit alpha
stiffvec = linspace(0, 4000, 25);
for j = 1: length(stiffvec)
    LD50fit0 = fit_params(2);
    p(1) = stiffvec(j);
    [ C, LD50fit, psimfit] = model_EMTpath( C_init, p, tsamp, LD50fit0, beta );
    LD50endfit(j) = LD50fit(end);
    
end
p(1) = stiff(3);
[ C, LD50mod, psimfit] = model_EMTpath( C_init, p, tsamp, LD50fit0, beta );
modelLD50 = LD50mod(end);
pct_error = 100*abs((modelLD50-outLD50(3))./(modelLD50))

figure;
subplot(1,2,1)
plot(stiff(1:2), outLD50(1:2), '*', 'LineWidth',2)
hold on 
plot(stiffvec, LD50endfit, 'LineWidth',2)
plot(stiff(3), outLD50(3), 'o', 'LineWidth', 5', 'color', 'r')
xlabel ('matrix stiffness \sigma (Pa)')
ylabel ('LD50 (\muM) after 6 days')
title ('A Model predicted versus measured LD50 MCF-7')
legend('measured LD50 for calibration', 'model R-sq = 0.74', 'measured 3rd stiffness % err = 46%')
legend boxoff
xlim ([0 stiffvec(end)])
ylim([ 0 22])
%%
outLD50(:,1) = [9.35, 12.09, 20.89];

params0 = [params(2) 8];% Initial Guess for alpha and LD50 base
options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
% test error function
%err = fit_model_EMT( params0,outLD50, stiff,C_init, const_params, tsamp, beta )
% Now fit to only two of the stiffnesses and LD50s
[fit_params, resnorm, reslsq]= lsqnonlin(@fit_model_EMT, params0, LB, UB, options, outLD50(1:2), stiff(1:2),C_init, const_params, tsamp, beta );
Rsq = 1- (sum((reslsq.^2))./(sum((mean(outLD50(1:2))-outLD50(1:2)).^2)))
% Model calibration
p = params;
p(2) = fit_params(1); % replace alpha with model fit alpha
stiffvec = linspace(0, 4000, 25);
for j = 1: length(stiffvec)
    LD50fit0 = fit_params(2);
    p(1) = stiffvec(j);
    [ C, LD50fit, psimfit] = model_EMTpath( C_init, p, tsamp, LD50fit0, beta );
    LD50endfit(j) = LD50fit(end);
    
end
p(1) = stiff(3);
[ C, LD50mod, psimfit] = model_EMTpath( C_init, p, tsamp, LD50fit0, beta );
modelLD50 = LD50mod(end);
pct_error = 100*abs((modelLD50-outLD50(3))./(modelLD50))

subplot(1,2,2)
plot(stiff(1:2), outLD50(1:2), '*', 'LineWidth',2)
hold on 
plot(stiffvec, LD50endfit, 'LineWidth',2)
plot(stiff(3), outLD50(3), 'o', 'LineWidth', 5', 'color', 'r')
xlabel ('matrix stiffness \sigma (Pa)')
ylabel ('LD50 (\muM) after 6 days')
title ('B Model predicted versus measured LD50 M6')
legend('measured LD50 for calibration', 'model R-sq =0.91 ', 'measured 3rd stiffness % error = 59% ')
legend boxoff
xlim ([0 stiffvec(end)])
ylim([ 0 22])



%% See how two parameters effect LD50 
% return to initial conditions from model above
C_i = C_init;
p  = params;

alphabase = C_init(
alphavec = linspace(0, , 5)
kdepvec =
for j = 1:length(sigmavec)
    p(1) = sigmavec(j);
for i = 1:length(TwG3vec)
    C_i(3) = TwG3vec(i);
    [ Cendtest, LD50eqtestmat(j,i), psimeqtest(j,i)] = model_EMTpath( C_i, p, tsamp, LD50base, beta );
end
end
%%
imagesc(LD50eqtestmat)
colorbar
ylabel ('mechanical stress')
xlabel('baseline concentration of TwG3')
title('Map of parameter space stress versus TwG3effect on LD50')

