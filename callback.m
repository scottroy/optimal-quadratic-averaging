function callback(x, f)
% adds f(x) to global variable callback_history (for plotting)

global callback_history;
callback_history(end+1) = f(x);

end

