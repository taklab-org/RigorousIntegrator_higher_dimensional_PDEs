clc
close all
clear

n = 4;%input('Select data: ');
switch n
    case 1 % First Equilibria
        load equilibria1_2D_SH.mat
        N  = [12,12];
        a = zeros(N(1)+1,N(2)+1);
        N_bar = size(a_bar);
        a(1:N_bar(1),1:N_bar(2)) = a_bar;
        L = [1,1.1];
        mu_k = mu_k_SH_sseig(N,3,L);
        mu_m = abs(mu_k(end,1));
        Nop_coeff = [0,0,0,-1];

        fprintf("\nFinding Equilibrium\n")
        [a,~] = newton_F_ss(a,L,mu_k,Nop_coeff);

        fprintf("\nFinding Eigenpair\n")
        DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff);
        
        [V,D] = eig(DF_ss); DD = flip(sort(diag(D)));
        [lambda,k] = max(diag(D));
        phi = reshape(V(:,k),N+1)/sum(V(:,k));  
        [phi,lambda,k] = newton_F_eig(a,phi,lambda,L,mu_k,Nop_coeff);


        fprintf("\nFinding both Equilibrium and Eigenpair\n")
        [a,phi,lambda,~] = newton_F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);
        


        %plot_cos_2d(a,phi,lambda,L);

    case 2
        load equilibria2_2D_SH.mat
        N  = [12,12];
        a = zeros(N(1)+1,N(2)+1);
        N_bar = size(a_bar);
        a(1:N_bar(1),1:N_bar(2)) = a_bar;
        L = [1,1.1];
        mu_k = mu_k_SH_sseig(N,3,L);
        mu_m = abs(mu_k(end,1));
        Nop_coeff = [0,0,0,-1];

        fprintf("\nFinding Equilibrium\n")
        [a,~] = newton_F_ss(a,L,mu_k,Nop_coeff);

        fprintf("\nFinding Eigenpair\n")
        DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff);
        
        [V,D] = eig(DF_ss); DD = flip(sort(diag(D)));
        [lambda,k] = max(diag(D));
        phi = reshape(V(:,k),N+1)/sum(V(:,k));  
        [phi,lambda,k] = newton_F_eig(a,phi,lambda,L,mu_k,Nop_coeff);


        fprintf("\nFinding both Equilibrium and Eigenpair\n")
        [a,phi,lambda,~] = newton_F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);

        plot_cos_2d(a,phi,lambda,L);

    case 3
        % zeros node is not a variable 

        load equilibria1_2D_DC.mat
        N  = [12,12];
        a = zeros(N(1)+1,N(2)+1);
        N_bar = size(a_bar);
        M1 = min(N_bar(1),N(1)+1);
        M2 = min(N_bar(2),N(2)+1);
        a(1:M1,1:M2) = a_bar(1:M1,1:M2) ;
        L = [1,1.1];
        epsilon = 0.4;
        sigma = 1;
        mu_k = mu_k_dcCH_sseig(N,epsilon,sigma,L);
        mu_m = abs(mu_k(end,1));
        Nop_coeff = [0,0,0,-1];

%         DF_fin = reshape(fin_dif_ss(a,L,mu_k,Nop_coeff),numel(a),numel(a));
%         DF = reshape(DF_steady_state(a,L,mu_k,Nop_coeff),numel(a),numel(a));
%         surf(abs(DF_fin - DF))
%         norm(DF_fin-DF)
%         



        

        fprintf("\nFinding Equilibrium\n")
        [a,~] = newton_F_ss(a,L,mu_k,Nop_coeff);
        if norm(a) < 10^-10
            error("Newton converged to the trivial solution.")

        end
        fprintf("\nFinding Eigenpair\n")
        DF_ss = reshape(DF_steady_state(a,L,mu_k,Nop_coeff),numel(a),numel(a));
        [V,D] = eig(DF_ss); DD = flip(sort(diag(D)));
        [lambda,k] = max(diag(D));
        phi = reshape(V(:,k),N+1)/sum(V(:,k));

        norm(F_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff))
        DF_eigenpair(a,phi,lambda,L,mu_k,Nop_coeff);

        [phi,lambda,k] = newton_F_eig(a,phi,lambda,L,mu_k,Nop_coeff);
        %
        %
                fprintf("\nFinding both Equilibrium and Eigenpair\n")
                [a,phi,lambda,~] = newton_F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);
%                 plot_cos_2d(a,phi,lambda,L);


    case 4
        load equilibria_3D_SH.mat
        N  = [4,4,4];
        a = zeros(N(1)+1,N(2)+1,N(3)+1);
        N_bar = size(a_bar);
        a(1:N_bar(1),1:N_bar(2),1:N_bar(3)) = a_bar;
        L = [1,1.1,1.2];
        mu_k = mu_k_SH_sseig(N,0.04,L);
        mu_m = abs(mu_k(end,1,1));
        Nop_coeff = [0,0,0,-1];



        fprintf("\nFinding Equilibrium\n")
        [a,~] = newton_F_ss(a,L,mu_k,Nop_coeff);

        fprintf("\nFinding Eigenpair\n")
        DF_ss = fin_dif_ss(a,L,mu_k,Nop_coeff);
        [V,D] = eig(DF_ss);  DD = flip(sort(diag(D)));
        [lambda,k] = max(diag(D));
        phi = reshape(V(:,k),N+1);
        [phi,lambda,k] = newton_F_eig(a,phi,lambda,L,mu_k,Nop_coeff);


        fprintf("\nFinding both Equilibrium and Eigenpair\n")
        [a,phi,lambda,~] = newton_F_ss_eig(a,phi,lambda,L,mu_k,Nop_coeff);
        %plot_cos_3d(a,phi,lambda,L)


end

nu = 1.3;
q = 2;
[Y0,Z0,Z1,Z2] = Bounds_steady_state(a,L,mu_k,Nop_coeff,nu,q,mu_m);
if Z1 >= 1
    error('Z_1 is greater than 1')
end
p = [Z2 ,-(1-Z0-Z1),Y0];
r = sort(roots(p));



DF_ss = DF_steady_state(a,L,mu_k,Nop_coeff);
[V,D] = eig(DF_ss); 
[DD,I] = sort(diag(D));
D_order = zeros(size(a));
P = zeros(numel(a),numel(a));
for i = 1:numel(a)
    P(:,i) = V(:,I(end-i+1))/(sum(abs(V(:,I(end-i+1)))));   
    D_order(i) = D(I(end-i+1),I(end-i+1));
end

[G_0,C,rho,epsilon,delta]  = Gershgorin_ring(a,P,mu_k,Nop_coeff,L,r(2),nu,mu_k(end,1));



























