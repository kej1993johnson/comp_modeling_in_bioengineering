% This code attempts to fit isolated growth and death rates for naive 231s
% after they respond to paclitaxel.

% We hope to use this type of experiment to test the effect of stiffness on
% sensitivity to drug by seeding cells in different matrix stiffnesses and
% 1. determining the growth and death rates (drug sensitivity) of the 
% population as a whole and also using isolated cell population
% trajectories (to fine ke km ge and gm) and then using these to estimate
% subpopulation composition and transition rates for each stiffness

% This data is of 231s dosed at 48 h for a total fo 72 hours. We will
% convert the percent confluence measurements to number of cells, and use
% the time after drug to find kpop and the time way after drug to find
% gpop. Then we will try fitting the data to the EMT model, assuming a low
% fraction of mesenchymal cells at the onset.
clear all; close all; clc
[N, T] = xlsread('../data/231_pac_example_data.xls');
time = N(:,1);
confluence = N(:,2);
cellnum = (confluence./100).*0.2e6;

%% First just graph the raw data

figure;
subplot(1,3,1)
plot(time, cellnum, '*', 'LineWidth',2)
xlabel ('time(hours)')
ylabel('Number of cells')
title('231 treatment 72h at 50 nM paclitaxel')

subplot(1,3,2)
plot(time, cellnum, '*', 'LineWidth',2)
xlabel ('time(hours)')
ylabel('Number of cells')
xlim([0 48])
title('231 pre-treatment')

subplot(1,3,3)
plot(time, cellnum, '*', 'LineWidth',2)
xlabel ('time(hours)')
ylabel('Number of cells')
xlim([120 445])
title('231 post-treatment')

%% Use data only from post treatment to regrowth before carrying capacity
ind1 = find(time >120);
ind2 = find(time <550);

t1= time(ind1);
t = t1(ind2);
n1= cellnum(ind1);
n = n1(ind2);

%% First find kill rate by fitting to single exponential death model
Nmin = min(n);
ind3 = find(n==Nmin);
td2 =t(1:ind3);
td= td2-td2(1); % make the first time analyzed t0
tg2 =t(ind3:end);
tg= tg2-tg2(1); % make the first time analyzed t0
nd = n(1:ind3);
ng = n(ind3:end);
N0g = ng(1);
N0d = nd(1);
exfitd = simmodeld(.005, td, N0d);
exfitg = simmodelg(.004, tg, N0g);

%% Bayesian parameter estimation to find growth and death rates separately
% single exponential
sigma = 0.1;
pfxform = @(pval)[1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1].*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

modelfund = @(p)simmodeld(p, td, N0d); % single exponential model with death
modelfung = @(p)simmodelg(p, tg, N0g); % single exponential model with death

dguess = -(yfxform(nd(end-5)) - yfxform(nd(5)))/(td(end-5)-td(5));
gguess = (yfxform(ng(end-5)) - yfxform(ng(5)))/(tg(end-5)-tg(5));

loglikelihoodd = @(phat)sum(log(normpdf(yfxform(nd),yfxform(modelfund(pbxform(phat))), sigma)));
loglikelihoodg = @(phat)sum(log(normpdf(yfxform(ng),yfxform(modelfung(pbxform(phat))), sigma)));
objfund = @(phat)-loglikelihoodd(phat);
objfung = @(phat)-loglikelihoodg(phat);
phatbestd = fminsearch(objfund, pfxform(dguess));
phatbestg = fminsearch(objfung, pfxform(gguess));
dfit = pbxform(phatbestd);
gfit = pbxform(phatbestg);

modelnd= simmodeld(dfit, td, N0d);
modelng= simmodelg(gfit, tg, N0g);

figure;
subplot(1,2,1)
plot(td, nd, '*', 'LineWidth', 2)
hold on
plot(td, modelnd, '-', 'LineWidth', 2)
xlabel('time (hours)')
ylabel('number of cells')
title(['231 drug response, d= ', num2str(dfit)])
legend ('data', 'model fit')
legend boxoff

subplot(1,2,2)
plot(tg, ng, '*', 'LineWidth', 2)
hold on
plot(tg, modelng, '-', 'LineWidth', 2)
xlabel('time (hours)')
ylabel('number of cells')
title(['231 drug response, g= ', num2str(gfit)])
legend ('data', 'model fit')
legend boxoff

% Make contcatenated data sets and models

figure;
plot(t, n, '*', 'LineWidth', 2)
xlabel('time(hours)')
ylabel('number of cells')
title('Example trajectory to be fit')

%% Now run some simulations of EMT model to see if you can generate a similar looking trajectory

% ODE model of EMT cells
% dE/dt = geE - kemE + kmeE -dE
% dM/dt = gmM - kmeM + kemE - dM

tsamp = t-t(1); % have time start at 0
N0 = N0d; % have initial cell number be right after dosing is removed
ge = 0.0057; % growth rate of epithelial cells
kem = 0.000032; % baseline transition rate E to M ** this will be a function of stiffness**
kme = 0.000025; % transition rate from M to E
gm = 0.0032; % growth rate of mesenchymal cells
fracE = 0.8;
de = 0.008;
dm= 1e-6;

p = [ ge; gm; kem; kme; de; dm; fracE];

C_init(1) = N0*fracE;
C_init(2) = N0*(1-fracE);

f = @(t,C) [ ge.*C(1) - kem.*C(1) + kme.*C(2)- de.*C(1); % dE/dt
             gm.*C(2) - kme.*C(2) + kem.*C(1)- dm.*C(2)]; % dM/dt

options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:2);
[t,C]=ode45(f, tsamp,C_init, options);
E=C(:,1);
M = C(:,2);
Nsim = E+M;

figure
subplot(1,2,1)
plot(tsamp, Nsim, '-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('number of cells')
title('Simulated E and M cell population response to drug')

subplot(1,2,2)
plot(tsamp, E, 'b-', 'LineWidth',2)
hold on
plot(tsamp, M, 'r-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('number of cells')
title('Simulated E and M cell population response to drug')
%% Now try fitting the above model to the dose response using
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);


% test function with parameters p

Nmeas=EMT_model(p, N0, tsamp) + 1000*rand(1, length(tsamp))';

%%
% kem, de, dm, fracE
p0 = [kem; de; dm; fracE];
pgiven = [ge; gm; kme];
%p = [ ge; gm; kem; kme; de; dm; fracE];
pguess = vertcat(pgiven(1:2), p0(1), pgiven(3), p0(2:4));

LB = vertcat(zeros*ones(length(p0)-1,1), 0);
UB = vertcat(Inf*ones(length(p0)-1,1),1);

[p_fit, resnorm, residuals] = lsqnonlin(@fit_EMT,...
    p0,...
    LB,...
    UB,...
    options,...
    pgiven,...
    tsamp,...
    Nmeas,...
    N0);
%% Find and evaluate model fit
 err_params = p_fit-p0;
 pfit_all = vertcat(pgiven(1:2), p_fit(1), pgiven(3), p_fit(2:4));
 
 Nmodel_fit=EMT_model(pfit_all, N0, tsamp);
 
 figure;
 plot(tsamp, Nmeas, '*')
 hold on 
 plot(tsamp, Nmodel_fit, '-', 'LineWidth',2)
 xlabel('time')
 ylabel('number of cells')
 title('Test fit to simulated data')
 %% Now perform fit to real data
 [p_fitdata, resnorm, residuals] = lsqnonlin(@fit_EMT,...
    p0,...
    LB,...
    UB,...
    options,...
    pgiven,...
    tsamp,...
    n,...
    N0);
 
pfit_data = vertcat(pgiven(1:2), p_fitdata(1), pgiven(3), p_fitdata(2:4));
 
 [Nmodel_fit_data, Efit, Mfit]=EMT_model(pfit_data, N0, tsamp);
 
 figure;
 subplot(1,2,1)
 plot(tsamp, n, '*')
 hold on 
 plot(tsamp, Nmodel_fit_data, '-', 'LineWidth',2)
 xlabel('time')
 ylabel('number of cells')
 title('Fit to real 231 data')
 subplot(1,2,2)
 plot(tsamp, n, '*')
 hold on 
 plot(tsamp, E, 'b-', 'LineWidth',2)
 plot(tsamp, M, 'r-', 'LineWidth', 2)
 xlabel('time')
 ylabel('number of cells')
 title('Fit to real 231 data E and M estimates')
 