% Simulate observed data X(t) = mu(t) + epsilon(t) where mu(t) = mu_0 +
% delta*H(t - t_0)

% Test hypothesis H0: delta = 0 vs H1: delta =/ 0 using likelihood ratio
% test
% set parameters
T = 50000; % length of time (time -> infinity)
t_0 = 25000; % time of change point corresponding to time
N = 50000; % number of simulations (N -> infinity)
mu_0 = 10; % mean 
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error

% Test statistic of likelihood ratio test asympotically follows a Chi
% square distribution with 1 degree of freedom 
critical_val = chi2inv(0.95,1); % 3.84 

% Repeat simulation with delta = 0 1000 times
delta = 0; % set delta
obs_delta_zero = cell(length(T), length(N)); % empty cell array to store test statistics 

% loop to generate data for each length of time with delta = 0 and sigma = 1
for i = 1:length(T)
    for j = 1:length(N)
        S = zeros(1, N(j)); % empty vector to store chi statistics 
        for k = 1:N(j)
            X = zeros(1, T(i)); % empty vector to store observed data 
            epsilon = normrnd(mu_epsilon, sigma, 1, T(i)); % generate random normal error 
            for t = 1:T(i)
                X(t) = mu_0 + delta*heaviside(t - t_0(i)) + epsilon(t);
            end
            var_0_hat = var(X, 1);
            var_1_hat = (1/T(i))*((t_0(i) - 1)*var(X(1:(t_0(i) - 1)), 1) + (T(i) - t_0(i) + 1)*var(X(t_0(i):T(i)), 1));
            % test statistic
            lambda_x = (var_1_hat/var_0_hat)^(T(i)/2);
            % -2*ln(lamda) approximately follows a Chi Square distribution with 1 df 
            chi_statistic = -2*log(lambda_x);
            S(k) = chi_statistic;
        end
        obs_delta_zero{i, j} = S;
    end
end

%% Check if false positive rate -> 0.05 

false_pos_rate = zeros(length(T), length(N));

% loop to calculate false positive rates 
for i = 1:size(false_pos_rate, 1)
    for j = 1:size(false_pos_rate, 2)
        false_pos_rate(i, j) = sum(obs_delta_zero{i, j} > critical_val)/length(obs_delta_zero{i, j});
    end
end

save('alpha_check.mat', 'obs_delta_zero', 'false_pos_rate')

exit 
