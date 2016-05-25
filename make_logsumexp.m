function [ objective ] = make_logsumexp(A, b, alpha)
% Makes log-sum-exp function
%   f(x) = log(sum(exp(Ax+b))) + 0.5*alpha*|x|^2

objective.f = @log_sum_exp;
objective.A = A;
objective.b = -b;
objective.gamma = alpha;

% no g

end

