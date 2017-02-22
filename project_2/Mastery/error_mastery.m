
N_vals = [16,32,64,128];% Values of N to solve over
[~,num] = size(N_vals);
E = zeros(1,num); % Array to store errors for each N value

for index = 1:num
    % Set N
    N = N_vals(index);
    % Construct x points
    dx = 2*pi/N;
    x = zeros(N,1);
    for j = 0:N-1
        x(j+1) = (j) * dx;
    end
    % Compute f(x)
    fx = sin(x);
    % Construct D2
    D2 = diag(2*ones(N,1)) + diag(-1*ones(N-1,1),1) + diag(-1*ones(N-1,1),-1);
    D2(N,1) = -1;
    D2(1,N) = -1;
    D2 = D2 / (dx*dx);
    % Construct P
    P = zeros(N);
    for i = 1:N
        if mod(i,2) == 1
            P((i + 1) / 2 , i) = 1;
        else
            P( N + 1 - i/2 ,i) = 1;
        end
    end
    PDP = P'*D2*P; % A symmetric positive semi-definite matrix
    Pfx = P'*fx;
    % Compute uN by solving (P^T D2 P) (P^T x) = P^T b
    PuN = algo_mastery(PDP,Pfx); 
    uN = P*PuN;
    % Compute Error
    E(index) = real(sqrt( dx * sum ( ( fx - uN ).^2 ) )); 
end

% Prints a table of errors
Error_Table = table(N_vals(:), E(:), 'VariableNames',{'N' 'Error'})
% Prints a table of Error multiplied by dx^2
scaledE = E.*(N_vals.^2);
ScaledError_Table = table(N_vals(:), scaledE(:), 'VariableNames',{'N' 'Scaled_Error'})