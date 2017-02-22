function [tp1 , tp4, tfft] = time_allalgos(N)
% Yadu Bhageria
% 00733164

% Takes as an input number of points N and returns the time taken by each
% of the three algos implemented

% Construct x points
dx = 2*pi/N;
x = zeros(N,1);
for j = 0:N-1
    x(j+1) = (j) * dx;
end
% Compute f(x)
fx = sin(x);
% Construct the A matrix for solving Au = f(x)
A = construct_A(N);

tic;
algo_p1(A,fx); % Time algo p1
tp1 = toc;

% Construct the Z matrix for solving - (Z Lambda Z') u = f(x)
Z = construct_Z(x,N);
% Construct the array storing the eigenvalues of D_2
lambda = construct_lambda(N);

tic;
algo_p4(Z,lambda,fx,N); % Time algo p4
tp4 = toc;

tic;
algo_fft(N,lambda,fx); % Time algo fft
tfft = toc;

end