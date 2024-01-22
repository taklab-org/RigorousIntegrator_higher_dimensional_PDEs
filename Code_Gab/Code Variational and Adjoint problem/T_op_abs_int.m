function [Ta] = T_op_abs_int(a)
N = size(a)-1;
Ta = intval(zeros(N+1));
otherdims = repmat({':'},1,length(N)-1);
for i = 1:N(end)-1
    Ta(otherdims{:}, i+1) =  abs(a(otherdims{:}, i+2)) + abs(a(otherdims{:}, i));
end
Ta(otherdims{:}, N(end)+1) =  abs(a(otherdims{:}, N(end)));
end

