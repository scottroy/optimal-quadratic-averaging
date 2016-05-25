function [a, b] = bracket_minimizer(f)
% Finds a < b such that [a,b] contains minimizer of one dimensional convex
% function f

% By convexity, it suffices to find a <= b with:
%       a <= aa for some aa satisfying f(a) >= f(aa)
%       b >= bb for some bb satisfying f(b) >= f(bb)

rho = 2;
aa = -1;
a = rho*aa;
faa = f(aa);
fa = f(a);
while faa > fa
    aa = a;
    faa = fa;
    
    a = rho*aa;
    fa = f(a);
end

bb = 1;
b = rho*bb;

fbb = f(bb);
fb = f(b);
while fbb > fb
    bb = b;
    fbb = fb;
    
    b = rho*bb;
    fb = f(b);
end


end

