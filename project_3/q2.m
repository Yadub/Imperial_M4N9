% Yadu Bhageria
% 00733164 

clear;
load('G.mat');  % Load given data
G = sparse(G);  % Make G sparse
[N,~] = size(G);% Set size of N
preconditioning = 0; % No preconditioning

b = ones(N,1);  % Initial column vector of ones

alphas = [0.5,0.7,0.9,0.99,0.9999];   % Values of alpha of interest
na = length(alphas);

x0 = zeros(N,1);% Initalize x0 for Power method
x0(1) = 1;      % Set x0 = e1

for i = 1:na
    alpha = alphas(i);          % Set alpha value
    A = speye(N) - alpha * G';  % Compute sparse matrix A
    tol = 1e-8;     % Set tolerance
    
    tic;
    [solg,count(i,1)] = GMRES(A, b, tol, preconditioning);   % GMRES
    tGMRES(i,1) = toc;  % Timetaken
    
    tic;
    [solp,countPM(i,1)] = sparsePageRank(G,alpha,tol,x0);     % Power method
    tPM(i,1) = toc;     % Timetaken
    
    
    solg = solg / norm(solg);
    solp = solp / norm(solp);
    norm(solg - solp)
    
    tol = 1e-10;     % Set new tolerance
    
    tic;
    [solg,count(i,2)] = GMRES(A, b, tol, preconditioning);   % GMRES
    tGMRES(i,2) = toc;  % Timetaken
    
    tic;
    [solp,countPM(i,2)] = sparsePageRank(G,alpha,tol,x0);     % Power method
    tPM(i,2) = toc;     % Timetaken
end

% Values of interest are stored in:
% count
% countPM
% format shorte;
% tGMRES
% tPM

% We can plot the values of the number of iterations required by each method
% and tolerance level. I found tables to be more useful
% fig_iters = figure();
% hold on;
% plot(alphas,count(:,1),'x');
% plot(alphas,countp(:,1),'o');
% plot(alphas,count(:,2),'x');
% plot(alphas,countp(:,2),'o');
% legend('1e-8 GMRES','1e-8 Power Method','1e-10 GMRES','1e-10 Power Method');
% hold off;
