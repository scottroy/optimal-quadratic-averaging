function [ val, grad ] = log_sum_exp(x)
% Returns value and gradient of function f(x) = log(sum(exp(x)))

 x_max = max(x);
 temp = exp(x - x_max);
 val = x_max + log(sum(temp));
 
 if nargout == 2
    grad = 1 / sum(temp) * temp;
 end

end

