close all, clear all, clc
[N1, T1] =xlsread('../data/YAP_data.xls');

% Make a structure to hold data
for i = 1:length(T1)
YAP(i).celltype = T1{i};
YAP(i).meanYAPlevel = N1(:,i);

ind = find(isnan(YAP(i).meanYAPlevel));
YAP(i).meanYAPlevel(ind)=[];
YAP(i).numcells = length(YAP(i).meanYAPlevel);
YAP(i).meantot = mean(YAP(i).meanYAPlevel);
YAP(i).CI = prctile(YAP(i).meanYAPlevel, [2.5 97.5]);
end

%% make figure to plot distributions
Color = {'b'; 'r';'g'; 'm';};

figure;
subplot(1,2,1)
for i = 1:2
histogram(YAP(i).meanYAPlevel,'BinWidth',2, 'FaceColor', Color{i})
hold on
end
xlabel ('Mean YAP nuclear level (A.U.s)')
ylabel('Frequency')
legend([YAP(1).celltype],[YAP(2).celltype] )
subplot(1,2,2)
for i = 3:4
histogram(YAP(i).meanYAPlevel,'BinWidth',100, 'FaceColor', Color{i})
hold on
legend([YAP(i).celltype])
end
xlabel ('Mean YAP nuclear level in each cell(A.U.s)')
ylabel('Frequency')
legend([YAP(3).celltype], [YAP(4).celltype])
%% Plot mean and CI
figure;
subplot(1,2,1)
for i = 1:2
plot(YAP(i).meantot, i, '*', 'LineWidth', 3,'color', Color{i})
hold on
plot(YAP(i).CI, [i,i],'-','LineWidth', 2,'color', Color{i})
hold on
end
ylim([0 4])
xlabel ('Mean YAP nuclear level all cells (A.U.s)')
title('Mean & 95 % CI MCF-7 hard and soft')
legend([YAP(1).celltype],[YAP(1).celltype],[YAP(2).celltype],[YAP(2).celltype] )
legend boxoff

subplot(1,2,2)
for i = 3:4
plot(YAP(i).meantot, i-2, '*', 'LineWidth', 3,'color', Color{i})
hold on
plot(YAP(i).CI, [i-2,i-2],'-','LineWidth', 2,'color', Color{i})
hold on
end
xlabel ('Mean YAP nuclear level (A.U.s)')
ylabel('Frequency')
title('Mean & 95 % CI 231s hard and soft')
ylim([0 4])
legend([YAP(3).celltype],[YAP(3).celltype], [YAP(4).celltype],[YAP(4).celltype])
legend boxoff
%% Next want to test if each distribution is statisticall significantly different from the other
% (within cell type I am assuming)
% can use a kstest to test if distributions are statistically significantly
% different... although I've heard this one is sort of an easy way to get
% to significance
for i = 1
[ h, p]= kstest2(YAP(i).meanYAPlevel, YAP(i+1).meanYAPlevel);
[h, pt] =ttest2(YAP(i).meanYAPlevel, YAP(i+1).meanYAPlevel);
YAP(i).pkstest = p;
YAP(i+1).pkstest = p;
YAP(i).pttest = pt;
YAP(i+1).pttest = pt;
end
for i = 3
[ h, p]= kstest2(YAP(i).meanYAPlevel, YAP(i+1).meanYAPlevel);
[h, pt] =ttest2(YAP(i).meanYAPlevel, YAP(i+1).meanYAPlevel);
YAP(i).pkstest = p;
YAP(i+1).pkstest = p;
YAP(i).pttest = pt;
YAP(i+1).pttest = pt;
end

