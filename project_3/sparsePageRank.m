function [ x, iteration_count ] = sparsePageRank( G, alpha, tol, x0 )
%PAGERANK Implementation of the PageRank Algorithm to find the PageRank
%state
%   Inputs: G: NxN  matrix of link transition probabilities from page i to j
%           alpha: probability of moving to a random page rather than
%           following a link on the current page
%           tol: the sensitivity of the algorithm. i.e. terminates when
%           |x^(k+1) - x^(k)| < tol
%           x0: Nx1 column vector containing the inital probability for
%           being on each webpage
%   Output: x: Nx1 column vector containing the PageRank state
%           iteration_count: number of iterations needed to attain x within the
%           tolerance

[N,~] = size(G); % Determine N

a = ~sum(G,2); % Initalize a. Old method from project 2 gave a small error

one = ones(N,1); % Store the one array to
PGT = alpha * sparse(G'); % Transpose of the G part of P

xold = zeros(N,1);  % Initalize xold
xnew = x0;          % For storing the xnew after each iteration
iteration_count = 0;% Set iteration count to 0

while norm(xnew - xold) > tol % Compute the PageRank state upto prescribed tolerance
    xold = xnew;    % Discard x old
    xnew = PGT * xold + ( alpha * one * (a' * xold) + (1-alpha) * one * (one' * xold) ) / N; % Compute x new
%     xnew = xnew / sum(xnew); % Normalize X new using the fact that sum(xnew) must equal 1 to make sure error is within eps
    iteration_count = iteration_count + 1; % Increment iteration count
end

x = xnew; % Set computed PageRank state as the output vector

end


