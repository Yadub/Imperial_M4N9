% Yadu Bhageria
% 00733164 

clear;
load('G.mat');  % Load given data
G = sparse(G);  % Make G sparse
[N,~] = size(G);% Set size of N
alpha = 0.85;   % Set alpha
tol = 1e-8;     % Set tolerance
preconditioning = 0; % No preconditioning

A = speye(N) - alpha * G';  % Compute sparse matrix A
b = ones(N,1);  % Initial column vector of ones

[sol,num_iters] = GMRES(A, b, tol, preconditioning);% Use GMRES to compute the solution
sol = sol / norm(sol); % Normalize the solution

[~,I] = sort(sol,'descend');
top50 = I(1:50);

num_iters % Outputs number of iterations taken

