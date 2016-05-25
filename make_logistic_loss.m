function [ objective ] = make_logistic_loss(X, y, alpha)
% Makes regularized logistic loss function
%   f(w) = (1 / N) * Sum_{i=1}^N log(1 + exp(-y_i x_i^T w)) + 0.5*alpha*|w|^2
% X is N x n data matrix (X(i,:) = x_i^T)
% y is N x 1 label vector (y(i) = y_i)
% alpha is regularization constant

objective.f = @mean_log_one_plus_exp;
objective.A = @(w) -y .* (X*w);
objective.AT = @(w) -X'*(y.*w);

objective.gamma = alpha;

% no b, no g
end

