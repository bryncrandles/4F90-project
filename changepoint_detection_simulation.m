% Simulate observed data X(t) = mu(t) + epsilon(t) where mu(t) = mu_0 +
% delta*H(t - t_0)

% Test hypothesis H0: delta = 0 vs H1: delta =/ 0 using likelihood ratio
% test
% set parameters
T = 200; % length of data
t_0 = 100; % time of change point
mu_0 = 10; % mean 
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error

% Test statistic of likelihood ratio test asympotically follows a Chi
% square distribution with 1 degree of freedom 
critical_val = chi2inv(0.95,1); % 3.84 

% Repeat simulation with delta = 0 1000 times
delta = 0; % set delta
obs_delta_zero = zeros(1, 1000); % empty vector to store test statistics 

% loop to generate data with delta = 0 and sigma = 1
for i = 1:1000
    X = zeros(1, T); % empty vector to store observed data 
    epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    var_hat_0 = var(X, 1);
    var_hat_1 = (1/T)*((t_0 - 1)*var(X(1:t_0 - 1), 1) + (T - t_0 + 1)*var(X(t_0:T)));
    % test statistic
    lambda_x = (var_hat_1/var_hat_0)^(T/2);
    % -2*ln(lamda) follows Chi Square distribution with (r0 - r) DF 
    chi_statistic = -2*log(lambda_x);
    obs_delta_zero(i) = chi_statistic;
end

false_pos_rate = sum(obs_delta_zero > critical_val)/length(obs_delta_zero);

% false positive rate is 0.025 = 2.5%

%% Repeat simulation 1000 times with delta = 0.01 to 10 using sigma = 1
delta = 0.01:0.01:5;  % set delta vector
obs_delta = zeros(length(delta), 1000); % empty matrix to store data

% loop to generate data with delta ranging from 0.01 to 10
for i = 1:length(delta)
    Y = zeros(1, 1000); % empty vector to store chi statistics 
    for j = 1:1000 
        X = zeros(1, T); % empty vector to store observed data 
        epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 
        for t = 1:T
            X(t) = mu_0 + delta(i)*heaviside(t - t_0) + epsilon(t);
        end
        var_hat_0 = var(X, 1);
        var_hat_1 = (1/T)*((t_0 - 1)*var(X(1:t_0 - 1), 1) + (T - t_0 + 1)*var(X(t_0:T), 1));
        % test statistic
        lambda_x = (var_hat_1/var_hat_0)^(T/2);
        chi_statistic = -2*log(lambda_x);
        Y(j) = chi_statistic;
    end
    obs_delta(i, :) = Y;
end

%% Plot power of the hypothesis test as function of delta
% (Recall power of a test: probability that test rejects H0 when H1 is true)
% power = 1 - false negative rate 

false_neg_rate_delta = sum(obs_delta < critical_val, 2)/length(obs_delta);
power = 1 - false_neg_rate_delta;

figure()
plot(delta, power)
title('Delta vs Power')
xlabel('Delta')
ylabel('Power')

%% Power of the hypothesis test as a function of delta and sigma 
sigma = 0:0.001:1;
delta = 0.01:0.1:5;
obs_sigma_delta = zeros(length(sigma), length(delta), 1000);

% loop to generate data for each combination of sigma and delta 
for i = 1:length(sigma)
    for j = 1:length(delta)
        Y = zeros(1, 1000);
        for k = 1:1000
            X = zeros(1, T); % empty vector to store observed data 
            epsilon = normrnd(mu_epsilon, sigma(i), T, 1); % generate random normal error 
            for t = 1:T
                X(t) = mu_0 + delta(j)*heaviside(t - t_0) + epsilon(t);
            end
            var_hat_0 = var(X, 1);
            var_hat_1 = (1/T)*((t_0 - 1)*var(X(1:t_0 - 1), 1) + (T - t_0 + 1)*var(X(t_0:T), 1));
            % test statistic
            lambda_x = (var_hat_1/var_hat_0)^(T/2);
            chi_statistic = -2*log(lambda_x);
            Y(k) = chi_statistic;
        end
        obs_sigma_delta(i, j, :) = Y;
    end
end

%% Plot power of test as function of sigma and delta 
false_neg_rate_sigma_delta = sum(obs_sigma_delta < critical_val, 3)/size(obs_sigma_delta, 3);
power = 1 - false_neg_rate_sigma_delta;



    