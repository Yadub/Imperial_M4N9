function [ P ] = construct_P( G, alpha)
%CONSTRUCT_P Constructs the matrix P from the matrix G and probability
%alpha
%   Inputs: G: NxN matrix containing link relations between webpages
%           alpha: probability of moving to a random page rather than
%           following a link on the current page
%   Output: PT: NxN matrix given in the assignment

[N,~] = size(G);% Determine N

% a = zeros(N,1); % Initalize a
% for i = 1:N
%     if sum(G(i,:)) < 1 - eps % Here machine epsilon is used for clarity
%         a(i) = 1; % Set values of a = 1 when the ith row of G doesn't sum to 1 (within error of machine epsilon)
%     end
% end
a = ~sum(G,2); % Does the same since rows of G sum either to 1 or are zero

one = ones(N,1); % Store the one array to
P = alpha * (G + a * one' / N ) + (1-alpha) * (one * one') / N; % Compute P

end

