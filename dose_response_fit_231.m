% This code laods in dose response data and fits them to the two population
% model

clear all, close all, clc

%pat = xlsread('../data/subpop_sorting_90_10.xls'); % replace 75 ADR 25 WT with 90:10
%pat = xlsread('../data/subpop_sorting_oldADR.xls');
data = xlsread('../data/dose_response_231.xls');
% make this into a structure
for i = 1:length(data)
    pat.stiff(i)= data(i,1);
    pat.dose(i)= data(i,2);
    pat.rep(i)=data(i,3);
    pat.viability(i)=data(i,4);
    pat.acc_time(i)=data(i,6);
   % pat.cell_type(i)= 'MCF-7';
end
%% Perform fitting
dose = data(:,2);
stiff = data(:,1);
rep = data(:,3);
viability = data(:,4);
acc_time= data(:,6);
isoft = find(stiff==200); % find indices of soft matrices
ihard = find(stiff ==2000); % find indices of hard matrices


%%
nsamp = 2;
params0 = horzcat( [25 0.01 25 0.01], 0.5.*ones(1, nsamp));%  LD50res, sloperes, LD50sens, slopesens, fres
paramslb = zeros( 1, 4+nsamp);
paramsub = horzcat( [ Inf 1 Inf 1], ones(1, nsamp));
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
               

    isoft = find(stiff==200); % find indices of soft matrices
    ihard = find(stiff ==2000); % find indices of hard matrices
    psingle0 = [ 25 0.1];

[Phard, resnormhard, residualshard] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    dose(ihard),...
    viability(ihard));
[Psoft, resnormsoft, residualssoft] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    dose(isoft),...
    viability(isoft));

[P, resnorm, residuals] = lsqnonlin(@fitmixedpops,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    viability,...
    nsamp);
%%
dmod = 0:1: max(dose);
model_soft = singlepopmodel( Psoft, dmod);
model_hard = singlepopmodel( Phard, dmod);
model_twopop(:,1)= twopopmodel2(P,P(5), dmod, nsamp);
model_twopop(:,2)= twopopmodel2(P,P(6), dmod, nsamp);
figure;
plot(dose(isoft), viability(isoft),'bo', 'LineWidth',2)
hold on
plot(dose(ihard), viability(ihard),'ro', 'LineWidth',2)
plot(dmod, model_soft, 'b-', 'LineWidth',2)
plot(dmod, model_hard, 'r-', 'LineWidth',2)
legend('200 Pa', '2000 Pa')
xlabel ('dose (\muM)')
ylabel('viability')
title('Single population model 231s')
legend boxoff


figure;
plot(dose(isoft), viability(isoft), 'bo', 'LineWidth',2)
hold on
plot(dose(ihard), viability(ihard), 'ro', 'LineWidth',2)
plot(dmod, model_twopop(:,1), 'b-', 'LineWidth',2)
plot(dmod, model_twopop(:,2), 'r-', 'LineWidth',2)
legend('200 Pa', '2000 Pa')
xlabel ('dose (\muM)')
ylabel('viability')
title('Two population model 231s')
legend boxoff


