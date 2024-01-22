clf

LW = 2;
% Fig 6(a)
load data_stripe/data_tau=0.00016923_OK2D_195steps.mat
% figure
semilogy(tau_vec(1:timestep-1),LineWidth=LW)
xlabel('$i (=1,2,\dots,K)$','interpreter','latex', 'FontSize', 30)
ylabel('Adaptive step size $\tau_i$','interpreter', 'latex', 'FontSize', 30)
set(gca,'FontSize',20)
% axis tight
% axis 'auto y'
% SaveFig(gcf,'Stepsize1')
% return

% Fig 6(b)
load data_spot/data_tau=3.2978e-05_OK2D_217steps.mat
figure
semilogy(tau_vec(1:timestep-1),LineWidth=LW,Color=[0.8500 0.3250 0.0980])
xlabel('$i (=1,2,\dots,K)$','interpreter','latex', 'FontSize', 30)
ylabel('Adaptive step size $\tau_i$','interpreter', 'latex', 'FontSize', 30)
set(gca,'FontSize',20)
% axis tight
% axis 'auto y'
% SaveFig(gcf,'Stepsize2')
% return

% Fig 6(c)
figure
load data_stripe/data_tau=0.00016923_OK2D_195steps.mat
semilogy(cumsum(tau_vec(1:timestep-1)),err_at_endpoint,LineWidth=LW), hold on
% 
load data_spot/data_tau=3.2978e-05_OK2D_217steps.mat
semilogy(cumsum(tau_vec(1:timestep-1)),err_at_endpoint,LineWidth=LW)
% 
xlabel('$t$','interpreter','latex', 'FontSize', 30)
ylabel('Error bounds \boldmath$\varrho_0$\boldmath','interpreter', 'latex', 'FontSize', 30)
yticks(fliplr([1 1e-2  1e-4 1e-6 1e-8]))
legend('Stripe pattern','Spot pattern','Location','southeast')
set(gca,'FontSize',20)
axis tight
hold off

% Stripe
% load data_stripe/data_tau=0.25_OK2D_30steps.mat
% semilogy(cumsum(tau_vec(1:timestep-1)),err_at_endpoint,LineWidth=3), hold on

% load data_stripe/data_tau=0.00011364_OK2D_136steps.mat
% semilogy(cumsum(tau_vec(1:timestep)),err_at_endpoint,LineWidth=3), hold on

% load data_stripe/data_tau=0.0002094_OK2D_143steps.mat
% semilogy(cumsum(tau_vec(1:timestep-1)),err_at_endpoint,LineWidth=3), hold on


% load data_stripe/data_tau=3.4686e-05_OK2D_310steps.mat
% semilogy(cumsum(tau_vec(1:timestep-1)),err_at_endpoint,LineWidth=3), hold on


% Spot
% load data_spot/data_tau=0.25_OK2D_36steps.mat
% semilogy(cumsum(tau_vec(1:timestep-3)),err_at_endpoint(1:end-2),LineWidth=3),hold on

% load data_spot/data_tau=0.00018946_OK2D_148steps.mat
% semilogy(cumsum(tau_vec(1:timestep)),err_at_endpoint,LineWidth=3)


% 