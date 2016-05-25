clear all; clc;
rng(0);

% make rosenbrock function
% n = 200;
% B = 10^6;
% alpha = 1;
% objective = make_rosenbrock(B, alpha);

% make logistic loss
X = load('a1a.txt');
y = X(:,1);
X = X(:,2:end);
[m, n] = size(X);
alpha = 10^(-4);
objective = make_logistic_loss(X,y,alpha);

% make log-sum-exp objective
% m = 100;
% n = 1000;
% A = rand(m,n);
% b = rand(m,1);
% alpha = 10^(-3);
% objective = make_logsumexp(A,b,alpha);
 
 
 
% set starting point and options
x0 = zeros(n,1);
options.tol = 10^(-3);
options.max_iters = 5000;
options.t = 1;
options.display = 'on';
options.linesearch_tol = 10^(-4);
options.display_interval = 25;


% Run and time
tic; result = optave(objective, alpha, x0, options); toc;
result