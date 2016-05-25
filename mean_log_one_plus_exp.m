function [ val, grad ] = mean_log_one_plus_exp(x)
% Returns value and gradient of function f(x) = mean(log(1 + exp(x)))

pos_ind = (x >= 0);

temp1 = 1 + exp(-x(pos_ind));
temp2 = 1 + exp(x(~pos_ind));

val =  mean([x(pos_ind) + log(temp1); log(temp2)]);

if nargout == 2
    grad = NaN(length(x),1);
    grad(~pos_ind) = 1 / length(x) * (1 - 1 ./ temp2);
    grad(pos_ind) = 1 / length(x) * 1 ./ temp1;
end

end

