function result = optave( objective, alpha, x, options )
% Implementation of the optimal quadratic averaging algorithm to minimize
% alpha-strongly convex function
%   F(x) = f(Ax - b) + g(x) + 0.5*gamma*norm(x)^2
% with beta-Lipschitz continuous gradient. Throughout F^* denotes the minimum
% value.
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
%   Giving F in a structured way (via a struct) is much better than just
%   providing a function handle for F.
%   The functions objective.A and objective.AT (assumed expensive) are
%   evaluated once per iteration, whereas objective.f and objective.g are
%   evaluated many times per iteration.
%   If you only provide a function handle for F, F will be evaluated many
%   times per iteration, which is usually very expensive.
%
%
%
% 
% alpha is either
%   1. strong convexity constant of F
%   2. *positive* lower bound on the strong convexity constant of F
%   
%
%
%
% x is the starting point for the algorithm
%
%
%
%
% options (not required) is a struct with fields
%
%
%  Field      Description                                   Default value
% ---------------------------------------------------------------------
%   t                   Memory size. Using memory (setting          10
%                       t >= 2) is only recommended in high
%                       dimensional problems, where the cost of
%                       solving a quadratic program in t
%                       variables is negligible.
%
%   tol                 Termination tolerance.                      10^(-3)
%
%   max_iters           Maximum iterations.                         1000
%
%   display             If options.display ='on', every             'on'
%   display_interval    display_interval iterations, the            25
%                       algorithm prints an upper bound on
%                       F^*, a lower bound on F^*, and the
%                       optimality gap.
%                       If t >= 2, the algorithm also prints
%                       the flag returned by the Matlab
%                       function quadprog (see Matlab
%                       documentation for what the flag
%                       means). If options.display = 'off',
%                       nothing is printed.
%
%   linesearch_tol      Termination tolerance for the linesearch.   10^(-4)
%                       The default setting works well, and you
%                       probably don't need to adjust this.
%                       If the algorithm stalls on a poorly
%                       conditioned problem, you might consider
%                       decreasing linesearch_tol.    
%
%   callback            callback is a function that is              None
%                       evaluated at xplus at the end of each
%                       iteration. Useful for collecting
%                       iteration data to plot.
%
%
%
%
% OUTPUT
%
%
% result is a struct with fields
%
%
%   Field              Description
% ---------------------------------------------------------------------
%   xplus_final        Final point xplus algorithm calculates. This is the
%                      approximate minimizer of F.
%
%   valxplus_final     Value of F at xplus_final. This is the approximate
%                      minimum value of F.
%
%   v_final            Final v algorithm calculates. This is a lower bound
%                      on the minimum value of F.
%
%   gap_final          valxplus_final - v_final. This is an upper bound on
%                      F(xplus_final) - F^*, where F^* is the minimum value
%                      of F.
%
%   num_iters          Number of iterations before the algorithm terminated.     
%
%
%
%


% make objective
[f, g, A, AT, b, gamma] = make_optave_objective(objective);

% set default options
if nargin == 3; options = []; end;
if ~isfield(options,'t'); options.t  = 10; end;
if ~isfield(options, 'tol'); options.tol = 10^(-3); end;
if ~isfield(options, 'max_iters'); options.max_iters = 1000; end;
if ~isfield(options, 'display'); options.display = 'on'; end;
if ~isfield(options, 'linesearch_tol'); options.linesearch_tol = 10^(-4); end;
if ~isfield(options, 'display_interval'); options.display_interval = 25; end;
callback = isfield(options, 'callback');

% set quadprog and fminbnd options
quadprog_options = optimoptions('quadprog', ... % used for optimal averaging
                                'Algorithm', 'active-set', ...
                                'Display', 'off');
                            
fminbnd_options = optimset('TolX', options.linesearch_tol); % used for line search


% evaluate F at x
Ax = A(x);
Ax_minus_b = b(Ax);
[valfx, gradfx] = f(Ax_minus_b);
[valgx, gradgx] = g(x);
valFx = valfx + valgx;
gradFx = AT(gradfx) + gradgx;
nx2 = norm(x)^2;
if gamma ~= 0
    valFx = valFx + 0.5*gamma*nx2;
    gradFx = gradFx + gamma*x;
end
AgradFx = A(gradFx);


% determine xplus by minimizing F on line between x and x - gradx    
dir = -gradFx;
ndir = norm(dir);
dir = dir / ndir;
Adir = -AgradFx / ndir;
ipxdir = x'*dir;
F_restricted = @(t) f(Ax_minus_b+t*Adir) + g(x + t*dir) + ...
    + 0.5*gamma*(nx2 + 2*t*ipxdir + t^2);
[t1, t2] = bracket_minimizer(F_restricted);
t_opt = fminbnd(F_restricted, t1, t2, fminbnd_options);
xplus = x + t_opt*dir;
nxplus2 = nx2 + 2*t_opt*ipxdir + t_opt^2;
Axplus = Ax + t_opt*Adir;
Axplus_minus_b = b(Axplus);
valxplus = f(Axplus_minus_b) + g(xplus) + 0.5*gamma*nxplus2;


c = x - gradFx / alpha;
Ac = Ax - AgradFx / alpha;
v = valFx - 0.5 * norm(gradFx)^2 / alpha;

gap = valxplus - v;

% set iteration counter
k = 0;

if options.t >= 2
% If options.t >= 2, store past t quadratic lower bounds
%   Q(x) = v_i + 0.5 * alpha * norm(x - x_i^{++})^2
% in vector V and matrix C: V(i) = v_i and C(:,i) = x_i^{++}
    V = v;
    C = c;
    AC = Ac; % AC stores vectors A*C(:,i)
    CTC = norm(c)^2; % CTC = C'*C
end

% print
if strcmpi(options.display, 'on')  
    % print options
    fprintf('Options:\n')
    fprintf('\tt = %d\n', options.t);
    fprintf('\ttol = %.3e\n', options.tol);
    fprintf('\tmax_iters = %d\n', options.max_iters);
    fprintf('\tlinesearch_tol = %.3e\n\n\n', options.linesearch_tol);
    
    if options.t >= 2
        heading_format = '%-5s | %-15s | %-15s | %-15s | %-5s \n';
        print_format = '%5d | %15.5f | %15.5f | %15.5f | %5s \n';
        
        fprintf(heading_format, 'Iter', 'upper', 'lower', 'gap', 'quadprog');
        fprintf(heading_format, '', 'bound', 'bound', '', 'flag');
    else
        heading_format = '%-5s | %-15s | %-15s | %-15s \n';
        print_format = '%5d | %15.5f | %15.5f | %15.5f \n';
        
        fprintf(heading_format, 'Iter', 'upper', 'lower', 'gap');
        fprintf(heading_format, '', 'bound', 'bound', '');
    end
    
    heading_separator = char(double('-')*ones(1,70));
    fprintf('%s\n', heading_separator);
    
    if options.t >= 2
        fprintf(print_format, k, valxplus, v, gap, '-');
    else
        fprintf(print_format, k, valxplus, v, gap);
    end
        
end

if callback; options.callback(xplus); end;

while gap > options.tol && k <= options.max_iters
    k = k + 1;
    
    % determine x by minimizing F on line between xplus and c
    dir = c - xplus;
    ndir = norm(dir);
    dir = dir / ndir;
    Adir = (Ac - Axplus) / ndir;
    ipxplusdir = xplus'*dir;
    F_restricted = @(t) f(Axplus_minus_b + t*Adir) + g(xplus + t*dir) + ...
        0.5*gamma*(nxplus2 + 2*t*ipxplusdir + t^2);
    [t1, t2] = bracket_minimizer(F_restricted);
    t_opt = fminbnd(F_restricted, t1, t2, fminbnd_options);
    x = xplus + t_opt*dir;
    nx2 = nxplus2+2*t_opt*ipxplusdir+t_opt^2;
    
    % evaluate F at x
    Ax = Axplus + t_opt * Adir;
    Ax_minus_b = b(Ax);
    [valfx, gradfx] = f(Ax_minus_b);
    [valgx, gradgx] = g(x);
    valFx = valfx + valgx;
    gradFx = AT(gradfx) + gradgx;
    if gamma ~= 0
        valFx = valFx + 0.5*gamma*nx2;
        gradFx = gradFx + gamma*x;
    end
    AgradFx = A(gradFx);
    
    
    % determine xplus by minimizing F on line between x and x-gradx
    dir = -gradFx;
    ndir = norm(dir);
    dir = dir / ndir;
    Adir = -AgradFx / ndir;
    ipxdir = x'*dir;
    F_restricted = @(t) f(Ax_minus_b+t*Adir) + g(x + t*dir) + ...
        + 0.5*gamma*(nx2 + 2*t*ipxdir + t^2);
    [t1, t2] = bracket_minimizer(F_restricted);
    t_opt = fminbnd(F_restricted, t1, t2, fminbnd_options);
    xplus = x + t_opt*dir;
    nxplus2 = nx2 + 2*t_opt*ipxdir + t_opt^2;
    Axplus = Ax + t_opt*Adir;
    Axplus_minus_b = b(Axplus);
    valxplus = f(Axplus_minus_b) + g(xplus) + 0.5*gamma*nxplus2;
    
    % calculate new quadratic lower bound
    vnew = valFx - 0.5 * norm(gradFx)^2 / alpha;
    cnew = x - gradFx / alpha;
    Acnew = Ax - AgradFx / alpha;
    
    % compute optimal averaging (lam,1-lam) of last two quadratics:
    %     Q(x) = vnew + 0.5 * alpha * norm(x - cnew)^2
    %     Q(x) = v + 0.5 * alpha * norm(x - c)^2
    d = norm(c-cnew)^2;
    if d == 0
        vmax = max(vnew, v);
        if vmax == vnew
            lam = 1;
        else
            lam = 0;
        end
    else
        lamHat = 0.5 + (vnew - v) / (alpha * d);
        lam = max(min(lamHat,1),0);
    end
    
    % if options.t == 1, calculate optimal averaging exactly
    if options.t == 1
        if lam == 1
            c = cnew;
            Ac = Acnew;
            v = vnew;
        else
            c = lam*cnew + (1-lam)*c;
            Ac = lam*Acnew + (1-lam)*Ac;
            v = v + 0.5*alpha*d*lam^2;
        end
        
    % if options.t >= 2, solve a quadratic program to determine
    % optimal averaging
    else
        % update CTC, C, AC, and V
        if k < options.t
            temp = C' * cnew;
            CTC = [CTC, temp; temp', norm(cnew)^2];
            C = [C, cnew];
            AC = [AC, Acnew];
            V = [V; vnew];
        else
            temp = C(:, 2:end)' * cnew;
            CTC = [CTC(2:end, 2:end), temp; temp', norm(cnew)^2];
            C = [C(:,2:end), cnew];
            AC = [AC(:,2:end), Acnew];
            V = [V(2:end); vnew];
        end
    
        % get dimension of averaging parameter lambda
        lambda_dim = length(V) + 1; 
        
        % set initial averaging parameter
        lambda0 = [zeros(lambda_dim-2,1);lam;1-lam];
        
        % compute optimal averaging
        temp = C' * c;
        H = alpha * [CTC, temp; temp', norm(c)^2];
        h = -0.5*diag(H)-[V;v];
        Aeq = ones(1,lambda_dim);
        beq = 1;
        lb = zeros(lambda_dim,1);
        ub = [];
        [lambda,v,qp_flag] = quadprog(H, h, [], [], Aeq, beq, lb, ub,...
                                      lambda0, quadprog_options);

        % update v, c, and Ac
        v = -v;
        c = [C,c]*lambda;
        Ac = [AC, Ac]*lambda;   
    end 
    
    % update gap
    gap = valxplus - v;
    
    % print
    if strcmpi(options.display, 'on') && mod(k,options.display_interval) == 0;
        if options.t >= 2
            fprintf(print_format, k, valxplus, v, gap, num2str(qp_flag));
        else
            fprintf(print_format, k, valxplus, v, gap);
        end
    end
    
    if callback; options.callback(xplus); end;

end

result.xplus_final = xplus;
result.valxplus_final = valxplus;
result.v_final = v;
result.gap_final = gap;
result.num_iters = k;

end

