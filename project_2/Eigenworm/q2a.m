% Yadu Bhageria
% CID: 00733164

format shorte;      % Set format as asked
load('theta.mat');  % Load given data

% Part (a)
C = theta' * theta; % Compute C

tic; % To see how long various implementations take
[W,T] = modified_house(C);
% [Q,T] = optimized_house(C);
timetaken = toc

d0 = diag(T);   % The leading diagonal of T
d1 = diag(T,1); % The first diagonal above (and below) of T
 
% part (b)
tol = 1e-14;                            % Set tolerance
lambda_w = QR_wilkinson(T, tol);        % Find eigenvalues with wilkinson shift
evals = sort(lambda_w, 'descend');      % Order lambda's by size

% It can be seen Rayleight shift gives very close values to Wilkinson shift from the code below
% lambda_r = QR_rayleighshift(T, tol);    % Find eigenvalues with rayleigh shift
% evals_r = sort(lambda_r, 'descend');    % Order lambda_r's by size
% lambda_diff = evals - evals_r;          % Compute difference
% lambda_error = norm(lambda_diff(1:10));     % Find error between the top 10 values


num_e = 10;         % Number of eigenvalues of interest
display('Most dominant eigenvalues (top10)');
display(evals(1:10)); % Prints 10 most dominant eigenvalues

% part (c)
Evecs = zeros(49,num_e);% Matrix to store the eigenvectors

for i = 1:num_e
    Evecs(:,i) = evec(C, evals(i)); % Compute the eigenvector associated with the ith largest eigenvalue
end

fig_evecs = figure(); % Plot the eigenvectors
hold on;
for i = 1:4
    plot(Evecs(:,i),'DisplayName',['Eval = ' num2str(evals(i))]);
end
legend(gca,'show');
title('Most dominant eigenvectors of C')
xlabel('Index');
ylabel('Value at Index');
hold off;

% To output the table
% pos = 1:49;
% data = [pos',Evecs(:,1),Evecs(:,2),Evecs(:,3),Evecs(:,4)]
% table.data = data
% table.dataFormat = {'%d', 1, '%.4e', 4}
% latexTable(table);

% part (d)
[N,~] = size(C);            % Set N
Theta_p = zeros(N,num_e);   % Initalize Theta_p
a = zeros(num_e,1);         % Initalize Array to hold the projection coefficients
Error = zeros(num_e,1);     % Initalize Array to hold the error
for p = 1:num_e
    a(p) = theta(1,:) * Evecs(:,p); % Find the projection coeffcient
    Theta_p(:,p) = Theta_p(:,max(p-1,1)) + a(p) * Evecs(:,p); % Project onto pth Principal Component
    Error(p) = sum( (theta(1,:) - Theta_p(:,p)').^2 ); % Compute the error
end

% Prints a table of projection coefficients
Proj_Table = table([1:num_e]', a(:), 'VariableNames',{'Principal_Component' 'Coefficient'})

fig_theta = figure;
hold on;
plot(theta(1,:),'--');
plot(Theta_p(:,2));
plot(Theta_p(:,4));
plot(Theta_p(:,6));
plot(Theta_p(:,10));
legend('Original Data','p = 2','p = 4','p = 6','p = 10');
title('Projected Data and Original Data');
hold off;

% Prints a table of errors
Error_Table = table([1:num_e]', Error(:), 'VariableNames',{'Principal_Component' 'Error'})
