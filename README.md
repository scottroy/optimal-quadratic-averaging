# optimal-quadratic-averaging
A Matlab implementation of the optimal quadratic averaging algorithm.
The accompanying paper can be found at http://arxiv.org/abs/1604.06543.

## Algorithm file descriptions:
The file ```optave.m``` contains a function that implements the optimal quadratic averaging algorithm.  A detailed description of ```optave.m``` is given in the comments in the function itself.

The function ```optave.m``` expects the objective in a structured form (see the comments in ```optave.m``` for more information).  The function in ```make_optave_objective.m``` is called by ```optave.m``` to process and check the structure of the objective passed into ```optave.m```. 

The file ```bracket_minimizer.m``` contains a method for bracketing the minimizer of a 1-dimensional convex function.  The function ```bracket_minimizer.m``` and the Matlab function ```fminbnd``` are used by ```optave.m``` to perform exact line search.

## Example file descriptions:
The file ```main.m``` shows examples of using ```optave.m``` to minimize the Rosenbrock function, the logistic loss function, and the log-sum-exp function.
The functions ```log_sum_exp.m```, ```mean_log_one_plus_exp.m```, and ```c_half_square_norm.m``` are component functions used to make the Rosenbrock function, the logistic loss function, and the log-sum-exp function.
The functions in ```make_rosenbrock.m```, ```make_logistic_loss.m```, and ```make_logsumexp.m``` put together the component functions and return structured objectives that can be passed into ```optave.m```.

Running the file ```memory_comparison_plot.m``` produces a plot that compares minimizing the logistic loss using ```optave.m``` for various memory sizes. 
The function ```callback.m``` is passed into ```optave.m``` to collect plotting data.

The files ```a1a.txt``` and ```colon-cancer.txt``` are LIBSVM data sets (for use in the logistic regression examples).
