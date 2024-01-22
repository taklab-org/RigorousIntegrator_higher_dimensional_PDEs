function cm = getting_cm(ba_sup,m,nu)
m1 = m(1); m2 = m(2); m3 = m(3);
%
psi = intval(zeros(m));
%
for k1=1:m1
  for k2=1:m2
    for k3=1:m3
      psi(k1,k2,k3) = psi_k(ba_sup,k1-1,k2-1,k3-1,m,nu);
    end
  end
end
%
cm = 3*wnorm(psi,nu);