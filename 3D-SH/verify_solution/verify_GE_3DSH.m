function success = verify_GE_3DSH(err_at_tau,ba_at_tau,nu1)
  success = 0;
  load trapping_region_3DSH_equilibria.mat a_mat C r rho nu
  a = a_mat;
  if nu ~= nu1
    error('Different nu')
  end
  if sup(err_at_tau)>(rho/C)
    error('Fail. Never proves the GE.')
  end
  N_bar = size(ba_at_tau);
  N = size(a);
  if any(N_bar>N)
    ia_ext = intval(zeros(N_bar));
    ia_ext(1:N(1),1:N(2),1:N(3)) = a;
    ba_ext = ba_at_tau;
  else
    ba_ext = zeros(N);
    ba_ext(1:N_bar(1),1:N_bar(2),1:N_bar(3)) = ba_at_tau;
    ia_ext = a;
  end
  a_tau_minus_ta = err_at_tau + wnorm(ba_ext-ia_ext,nu1) + r(2);
  disp(['Distance from target equiliburium: ',num2str(sup(a_tau_minus_ta))])
  if a_tau_minus_ta < (rho/C)
    success = 1;
  end
end