% Yadu Bhageria
% 00733164

N_vals = [16,32,64,128];% Values of N to solve over
[~,num] = size(N_vals);
E = zeros(1,num); % Array to store errors for each N value

for index = 1:num
    
    N = N_vals(index);
    
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
    % Solve for u using QR factorization (householder triangularization)
    u = algo_p1(A,fx);
    
    % hold on; plot(x,u); % To visually test the output
    
    E(index) = sqrt( dx * sum ( ( fx - u ).^2 ) );
end

% Prints a table of errors
Error_Table = table(N_vals(:), E(:), 'VariableNames',{'N' 'Error'})
% Prints a table of Error multiplied by dx^2
scaledE = E.*(N_vals.^2);
ScaledError_Table = table(N_vals(:), scaledE(:), 'VariableNames',{'N' 'Scaled_Error'})