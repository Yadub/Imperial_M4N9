% Yadu Bhageria
% 00733164 

clear;
load('G.mat');  % Load given data
G = sparse(G);  % Make G sparse
[N,~] = size(G);% Set size of N
preconditioning = 1; % With preconditioning

b = ones(N,1);  % Initial column vector of ones

alphas = [0.5,0.7,0.9,0.99,0.9999,0.999999];   % Values of alpha of interest
na = length(alphas);

x0 = zeros(N,1);% Initalize x0 for Power method
x0(1) = 1;      % Set x0 = e1

for i = 1:na
    alpha = alphas(i);          % Set alpha value
    A = speye(N) - alpha * G';  % Compute sparse matrix A
    tol = 1e-8;     % Set tolerance
    
    [solnp1,countnp(i,1)] = GMRES(A, b, tol, 0);   % GMRES without preconditioning
    [solp1,count(i,1)] = GMRES(A, b, tol, preconditioning);   % GMRES

    solnp1 = solnp1 / norm(solnp1);
    solp1 = solp1 / norm(solp1);
    diffe8(i) = norm(solp1 - solnp1);
    
    tol = 1e-10;     % Set new tolerance
    [solnp,countnp(i,2)] = GMRES(A, b, tol, 0);   % GMRES without preconditioning
    [solp,count(i,2)] = GMRES(A, b, tol, preconditioning);   % GMRES
    
    solnp = solnp / norm(solnp);
    solp = solp / norm(solp);
    diffe10(i) = norm(solp - solnp);
end

% Values of iterest
% count
% countnp
% format shorte;
% diffe8
% diffe10