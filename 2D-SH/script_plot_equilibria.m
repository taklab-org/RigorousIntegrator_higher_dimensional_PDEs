% Figure 2(a)
clear
load('verify_solution/trapping_region_2DSH_equilibria1.mat')
addpath('defect/')
L1 = 1; L2 = 1.1;
figure
plot_SH_profile(mid(a_mat),[L1,L2])
view([0 90])
colorbar
% SaveFig(gcf,'2DSH_equilibria1')

% Figure 2(b)
load('verify_solution/trapping_region_2DSH_equilibria2.mat')
figure
plot_SH_profile(mid(a_temp),[L1,L2])
view([0 90])
colorbar
% SaveFig(gcf,'2DSH_equilibria2')


% Figure 3(a)
clear
load data_stripe/a_bar_SH2D_50steps.mat
L1 = 1; L2 = 1.1;
figure
plot_SH_profile(mid(a_end),[L1,L2])
view([0 90])
colorbar
% SaveFig(gcf,'Init_2DSH_equilibria1')

% Figure 3(b)
load('verify_solution/trapping_region_2DSH_equilibria1.mat','a_mat')
a_tmp = zeros(size(a_end)); a_tmp(1:size(a_mat,1),1:size(a_mat,2))=mid(a_mat);
figure
plot_SH_profile(a_tmp - a_end,[L1,L2])
view([0 90])
colorbar
% SaveFig(gcf,'diff_2DSH_equilibria1')


% Figure 4(a)
clear
load data_spot/a_bar_SH2D_300steps.mat a_end
L1 = 1; L2 = 1.1;
figure
plot_SH_profile(mid(a_end),[L1,L2])
view([0 90])
colorbar
% SaveFig(gcf,'Init_2DSH_equilibria2')

% Figure 4(b)
load('verify_solution/trapping_region_2DSH_equilibria2.mat','a_temp')
a_tmp = zeros(size(a_temp)); a_tmp(1:size(a_end,1),1:size(a_end,2))=(a_end);
figure
plot_SH_profile(a_temp - a_tmp,[L1,L2])
view([0 90])
colorbar
% SaveFig(gcf,'diff_2DSH_equilibria2')
