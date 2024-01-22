function [n] = norm_1_nu_int(a,nu)
N = size(a)-1;

nu_vect = intval(zeros(N+1));
% otherdims = repmat({':'},1,length(N)-1);
% for i = 0:N(end)
%     if i ==0 
%     nu_vect1(otherdims{:}, i+1) =  1;
%     else
%     nu_vect1(otherdims{:}, i+1) =  nu^i;
%     end
% end

for i = 1:numel(mid(a))
    k = iso_coeff(i,N);
    nu_vect(i) = prod(2.^double(k>0))*nu^(k(end));
end



n = abs(a).*nu_vect;

for i = 1:length(N)
    n = sum(n);
end
end