function [u] = algo_p4(Z,lambda,fx,N)
% Yadu Bhageria
% 00733164

% Solves the linear system - Z Lambda Z^-1 u = fx with a periodic f and
% while enforcing the zero-mean condition
%   Inputs: matrix Z with the eigenvectors of D_2
%           array lambda with the eigenvalues of D_2
%           array fx with values of the periodic zero-mean function f(x)
%   Outputs: array u satisfying the equation - D_2 u = fx with the
%            zero-mean and periodic conditions


% -> u = - Z inv(Lambda) Z' f(x)
u = Z'*fx;
% Now multiply by inv(Lambda)
u(2:N) = u(2:N) ./ lambda(2:N); % divide by eigenvalue

% Finally multiply by Z and take the negative to obtain the solution
u = - Z*u;

end
