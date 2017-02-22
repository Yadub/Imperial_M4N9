function [u] = algo_fft(N,lambda,fx)
% Yadu Bhageria
% 00733164

% Solves -D_2*u = f using the fast fourier transform method with perodic
% boundary conditions and mean-value condition.
% Input: N: integer stating size of the problem
%        lambda: eigenvalues of D_2
%        fx: Nx1 values of the function f(x)
% Outputs: u: mx1 column vector containing the solution

% -> u = - Z inv(Lambda) Z' f(x)
u = fft(fx);
% Now multiply by inv(Lambda)
u(2:N) = u(2:N) ./ lambda(2:N); % divide by eigenvalue

% Finally multiply by Z and take the negative to obtain the solution
u = - ifft(u);

end