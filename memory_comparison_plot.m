clear all; clc;

% load data set and make logistic loss
X = load('a1a.txt');
y = X(:,1);
X = X(:,2:end);
[m, n] = size(X);
alpha = 10^(-4);
objective = make_logistic_loss(X,y,alpha);
plot_title = 'Minimizing Logistic Loss on LIBSVM a1a Dataset with $\alpha = 10^{-4}$';


% set starting point and options
x0 = zeros(n,1);
options.tol = 10^(-9);
options.max_iters = 5000;


% Run optave for multiple values of t and plot
global callback_history;
[f, g, A, AT, b, gamma] = make_optave_objective(objective);    
F = @(x) f(b(A(x))) + g(x) + 0.5*gamma*norm(x)^2;
options.callback = @(x) callback(x,F);
options.display = 'off';
    
xplus_values = {};
MIN = Inf;
mem_sizes = [1,2,5,10,15];
line_colors = 'rbgkm';
for i = 1:length(mem_sizes)
    fprintf('Running optave with t = %d.\n', mem_sizes(i));
    
    callback_history = [];
    options.t = mem_sizes(i);
    tic; result = optave(objective,alpha,x0,options); toc;
    
    xplus_values{end+1} = callback_history;
    MIN = min(MIN, min(callback_history));
end

% plot data
for i = 1:length(mem_sizes)
    semilogy(1:length(xplus_values{i}), xplus_values{i} - MIN, line_colors(i));
    hold on;
end
hold off;
title(plot_title, 'interpreter', 'latex');
legend({'$t=1$', '$t=2$', '$t=5$', '$t=10$', '$t=15$'},'interpreter', 'latex');
ylabel('$F(x_k^+) - F^*$', 'interpreter', 'latex');
xlabel('$k$','interpreter', 'latex');