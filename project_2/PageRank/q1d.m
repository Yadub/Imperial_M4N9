% Yadu Bhageria
% CID: 00733164

format shorte;  % Set format
load('G.mat');  % Load given data
alpha = 0.85;   % Set alpha
tol = 1e-8;     % Set tolerance

[N,~] = size(G);% Determine N

x0 = zeros(N,1);% Initalize x0
x0(1) = 1;      % Set x0 = e1

% P = construct_P(G, alpha); % Compute the transition matrix P
tic;
[x, count] = sparsePageRank(G,alpha,tol,x0); % Perform the PageRank algorithm
timetaken = toc;
[~,I] = sort(x,'descend'); % Compute which pages i correspond to the largest probabilties
top50 = I(1:50); % Store the indicies of the top 50 pages with largest probabilities

timetaken, count, top50(1:3)