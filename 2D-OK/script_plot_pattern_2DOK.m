% Fig 5(a)
load('../Code_Gab/Code Steady State and Eigenvalue Problem Specific Fm/equilibria1_2D_DC.mat')
% 
L1 = 1; L2 = 1.1; L = [L1,L2];
addpath('defect/')
figure
plot_DC_pattern(a_bar,L)
view([0 90])
% SaveFig(gcf,'2DOK_pattern1')

% Fig 5(b)
load('../Code_Gab/Code Steady State and Eigenvalue Problem Specific Fm/equilibria2_2D_DC.mat')
% 
L1 = 1; L2 = 1.1; L = [L1,L2];
% addpath('defect/')
figure
plot_DC_pattern(a_bar,L)
view([0 90])
% SaveFig(gcf,'2DOK_pattern2')
