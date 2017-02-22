function [ x, iteration_count ] = PageRank( G, alpha, tol, x0 )
%PAGERANK Implementation of the PageRank Algorithm to find the PageRank
%state
%   Inputs: G: NxN matrix of link transition probabilities from i to j
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

P = construct_P( G, alpha); % Construct P

PT = P'; % Store the transpose of P

xold = zeros(N,1);  % Initalize xold
xnew = x0;          % For storing the xnew after each iteration
iteration_count = 0;% Set iteration count to 0

while norm(xnew - xold) > tol % Compute the PageRank state upto prescribed tolerance
    xold = xnew;    % Discard x old
    xnew = PT * xold; % Compute x new
    xnew = xnew / sum(xnew); % Normalize X new using the fact that sum(xnew) must equal 1
    iteration_count = iteration_count + 1; % Increment iteration count
end

x = xnew; % Set computed PageRank state as the output vector

end

