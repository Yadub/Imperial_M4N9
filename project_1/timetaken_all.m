% Yadu Bhageria
% 00733164

% Compute array of N values to iterate over
max_power = 10; 
x = 4:max_power;
N_vals = 2.^x;
num = length(N_vals);

% Initalize arrays to store time taken for each method
tp1 = zeros(num,1);
tp4 = zeros(num,1);
tfft = zeros(num,1);

for i = 1:num
    N = N_vals(i); % Set N
    [tp1(i), tp4(i), tfft(i)] = time_allalgos(N); % Compute time taken
end

loglog(N_vals,tp1,N_vals,tp4,N_vals,tfft); % Plot on a loglog scale
legend('Algo 1', 'Algo 4', 'Algo fft'); 
title('Yadu Bhageria: Time taken comaprison between algo p1, algo p4, and algo fft');