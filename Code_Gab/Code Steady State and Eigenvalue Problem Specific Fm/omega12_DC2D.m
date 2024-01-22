function [N1,N2] = omega12_DC2D(N_Fm,N_Fm_ext,nu,q)

a_ones = intval(ones( sum(N_Fm_ext+1) ,1));
a_ones = iso_vec2mat_Fm_int(a_ones,N_Fm_ext);

Ma = size(a_ones)-1;
[m,n] = meshgrid(0:Ma(2),0:Ma(1));
alpha = 4*intval(ones(size(a_ones)));
alpha(1,:) = 2;
alpha(:,1) = 2;
alpha(1,1) = 1;

omega_k = (1+n+m).^q.*alpha.*nu.^(m+n);
omega_k = iso_mat2vec_Fm_int(omega_k,N_Fm_ext);

DF = omega_k;


N1 = [];
N2 = [];



N_Fm_pad = [N_Fm,-ones(1,length(N_Fm_ext) - length(N_Fm)) ];


k = 1;
i_ind = 0;
for i = N_Fm_ext
    for n = 0:i
        if n <= N_Fm_pad(i_ind+1) 
            N1 = [N1,DF(k)];
        elseif n > N_Fm_pad(i_ind+1) 
            N2 = [N2,DF(k)];
        end
        k = k+1;
    end
    i_ind =  i_ind  +1 ;
end

end

