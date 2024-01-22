function [Lambda] = lambda_op_int(a)
N = size(a)-1;
Lambda = intval(zeros(N+1));
otherdims = repmat({':'},1,length(N)-1);
for i = 0:N(end)
    Lambda(otherdims{:}, i+1) =  2*(i)*a(otherdims{:}, i+1) ;
end
end

