function [norm_B] = norm_B_nu_int(B,nu,N)
B = reshape(abs(B),[sqrt(numel(inf(B))),sqrt(numel(inf(B)))]);
M = size(B);
vect = intval(zeros(M(1),1));
for i = 1:M(1)
vect(i) = norm_1_nu_int(reshape(B(i,:),N+1),nu);
end

nu_vect = intval(zeros(N+1));
otherdims = repmat({':'},1,length(N)-1);
for i = 0:N(end)
    nu_vect(otherdims{:}, i+1) =  nu^(-i);
end

norm_B = reshape(vect,N+1).*nu_vect;

for i = 1:length(N)
    norm_B = max(norm_B);
end

end

