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
    for j = 1:N
        x(j) = (j-1) * dx;
    end

    % Compute f(x)
    fx = sin(x);
    % Construct the Z matrix for solving - (Z Lambda Z') u = f(x)
    Z = construct_Z(x,N);
    % Construct the array storing the eigenvalues of D_2
    lambda = construct_lambda(N);
    % Solve the system using the algorithm from part 4
    u = algo_p4(Z,lambda,fx,N);

    % hold on; plot(x,u); % To visually test the output
    
    E(index) = sqrt( dx * sum ( ( fx - u ).^2 ) ); % Store the error value
end

% Prints a table of errors
Error_Table = table(N_vals(:), E(:), real(E(:)), 'VariableNames',{'N' 'Error', 'Real_Error'})
% Prints a table of Error multiplied by dx^2
scaledE = E.*(N_vals.^2);
ScaledError_Table = table(N_vals(:), scaledE(:), real(scaledE(:)), 'VariableNames',{'N' 'Scaled_Error', 'Real_Scaled_Error'})