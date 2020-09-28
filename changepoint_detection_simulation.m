% Simulate observed data X(t) = mu(t) + epsilon(t) where mu(t) = mu_0 +
% delta*H(t - t_0)

% set parameters
T = 200; % length of data
t_0 = 100; % time of change point
mu_0 = 10; % mean 
mu_epsilon = 0; % mean of error
sigma = 0.1; % variance of error

% Repeat simulation with delta = 0 100 times
delta = 0; % set delta
y_0 = zeros(1, 100); % empty vector to store test statistics 

% loop to generate data with delta = 0
for i = 1:100 
    X = zeros(1, T); % empty vector to store obs data 
    epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    % Test hypothesis H0: delta = 0 vs H1: delta =/ 0 using likelihood ratio
    % test
    mu_hat_0 = mean(X(1:t_0 - 1));
    mu_hat_1 = mean(X(t_0:T));
    delta_hat = mu_hat_1 - mu_hat_0;
    var_hat_0 = var(X);
    var_hat_1 = var(X(1:t_0 - 1)) + var(X(t_0:T));
    % test statistic
    lambda_x = (var_hat_1/var_hat_0)^(T/2);
    % -2*ln(lamda) follows Chi Square distribution with (r0 - r) DF 
    chi_statistic = -2*log(lambda_x);
    y_0(i) = chi_statistic;
end

% Repeat simulation 100 times with delta = 10
delta = 10;  % set delta 
y_10 = zeros(1, 100);

% loop to generate data with delta = 10
for i = 1:100 
    X = zeros(1, T); % empty vector to store obs data 
    epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    % Test hypothesis H0: delta = 0 vs H1: delta =/ 0 using likelihood ratio
    % test
    mu_hat_0 = mean(X(1:t_0 - 1));
    mu_hat_1 = mean(X(t_0:T));
    delta_hat = mu_hat_1 - mu_hat_0;
    var_hat_0 = var(X);
    var_hat_1 = var(X(1:t_0 - 1)) + var(X(t_0:T));
    % test statistic
    lambda_x = (var_hat_1/var_hat_0)^(T/2);
    % -2*ln(lamda) follows Chi Square distribution with (r0 - r) DF 
    chi_statistic = -2*log(lambda_x);
    y_10(i) = chi_statistic;
end

% critical value
critical_val = chi2inv(0.95,1); % 3.84 

% When delta = 0, do not reject H0: delta = 0 at alpha = 0.05 for all
% simulations. When delta = 10, reject H0: delta = 0 at alpha = 0.05 for all 
% simulations.