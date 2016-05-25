function [ objective ] = make_rosenbrock(B, alpha)
% Makes Rosenbrock function on R^n
%   f(x) = 0.5*B*( (1-x(1))^2 + + Sum_{i=1}^{n-1} (x(i) - x(i+1))^2 + x(n)^2 ) 
%           + 0.5*alpha*|x|^2

objective.f = @(x) c_half_square_norm(x, B);
objective.A = @(x) [-x(1); -diff(x); x(end)];
objective.AT = @(x) diff(x);
objective.b = @(x) [x(1)+1;x(2:end)];

objective.gamma = alpha;

% no g

end

