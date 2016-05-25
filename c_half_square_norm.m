function [ val, grad ] = c_half_square_norm(x, c)
% Returns value and gradient of function f(x) = 0.5*c*norm(x)^2

val = 0.5*c*norm(x)^2;
if nargout == 2
    grad = c*x;
end

end

