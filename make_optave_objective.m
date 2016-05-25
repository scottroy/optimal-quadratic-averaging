function [f, g, A, AT, b, gamma] = make_optave_objective(objective)
% Makes objective
%   F(x) = f(Ax - b) + g(x) + 0.5*gamma*norm(x)^2
% for the optave algorithm.
%
% Checks objective and sets defaults for f, g, A, AT, b, and gamma.
% Converts matrix A to function handle for linear operator.
% Converts vector b to function handle for shift operator.
%
%
% INPUT
%
%
% objective is either
% 
%   1. function handle for F (to minimize)
%   2. struct of f, A, b, g, and gamma where F (to minimize) is
%           F(x) = f(Ax - b) + g(x) + 0.5*gamma*norm(x)^2
%      More precisely,
%      - objective.f is function handle for f
%        (f is required, no default)
%      - objective.A is a matrix or function handle for a linear operator
%        If objective.A is a function handle, you must also give
%        objective.AT, the function handle for the adjoint of A
%        (By default, A = I)
%      - objective.b is a vector or a function handle for a shift
%        operator x -> x - b
%        (By default, b = 0)
%      - objective.g is function handle for g
%        (By default, g = 0)
%      - objective.gamma is nonnegative number gamma
%        (By default, gamma = 0)
%
%
% OUTPUT
%
% Function handles for f, g, A (linear operator), AT (adjoint of A), b
% (shift operator), and number gamma
%



% define default functions
function [gx, gradgx] = default_g(x)
    gx = 0;
    gradgx = 0;
end

default_A = @(x) x;
default_AT = @(x) x;
default_b = @(x) x;


% if objective is a function handle, make it into a struct
if isa(objective, 'function_handle')
    temp = objective;
    clear objective;
    objective.f = temp;
end

% set f
if ~isfield(objective, 'f')
    error('You must provide the function handle objective.f.');
else
    f = objective.f;
end

% set g
if isfield(objective, 'g')
    g = objective.g;
else
    g = @default_g;
end

% set A and AT
if isfield(objective, 'A')
    if isa(objective.A, 'function_handle')
        A = objective.A;
        if ~isfield(objective, 'AT')
            error('You must provide a function handle for the adjoint of objective.A in objective.AT.');
        else
            AT = objective.AT;
        end
    else % A is a matrix
        A = @(x) objective.A*x;
        AT = @(x) objective.A'*x;
    end
else
    A = default_A;
	AT = default_AT;
end

    
% set b (b is shift operator x -> x - b)
if isfield(objective, 'b')
    if isa(objective.b, 'function_handle')
        b = objective.b;
    else % b is vector
        b = @(x) x - objective.b;
    end
else
    b = default_b;
end

% set gamma
if ~isfield(objective, 'gamma')
    gamma = 0;
else
    gamma = objective.gamma;
end
    

end

