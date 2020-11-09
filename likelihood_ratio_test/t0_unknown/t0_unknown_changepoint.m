%% Consider the case when t_0 is considered to be unknown.

% First try T = 200, code in t_0 = 100, delta = 10. 
% Test to see if the likelihood function is maximized at 
% t_0 = 100.

T = 200; % length of time series
t_0 = 100; % time of change point
mu_0 = 10; % mean of time series 
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error
epsilon = normrnd(mu_epsilon, sigma, 1, T); % error vector
delta = 5; % magnitude of change in the mean 
N = 1000;
critical_val = chi2inv(0.95,1); % 3.84 

%% check if t0 is detected when delta = 5
S = zeros(1, N); % empty vector to store test statistics 
t_0_hat_val = zeros(1, N); % empty vector to store t_0_hat values 

%% loop to generate time series of length T
for i = 1:N
    X = zeros(1, T); % initialize vector 
    epsilon = normrnd(mu_epsilon, sigma, 1, T); % error vector
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    var_0_hat = var(X, 1);
    t_0_val = 1:200; % vector of t_0 values to test
    % initialize vectors to store values of mu_0_hat, delta_hat, var_hat
    % and L (likelihood) that correspond to each t_0 value
    mu_0_hat_val = zeros(1, length(t_0_val));
    delta_hat_val = zeros(1, length(t_0_val));
    var_1_hat_val = zeros(1, length(t_0_val));
    L_val = zeros(1, length(t_0_val));
    for j = 1:length(t_0_val)
        mu_0_hat = mean(X(1:(t_0_val(j) - 1)));
        mu_0_hat_val(j) = mu_0_hat;
        xbar_1_hat = mean(X(1:(t_0_val(j) - 1)));
        xbar_2_hat = mean(X(t_0_val(j):T));
        delta_hat = xbar_2_hat - xbar_1_hat;
        delta_hat_val(j) = delta_hat;
        var_1_hat = (1/T)*((t_0_val(j) - 1)*var(X(1:(t_0_val(j) - 1)), 1) + (T - t_0_val(j) + 1)*var(X(t_0_val(j):T), 1));
        var_1_hat_val(j) = var_1_hat;
        L_val(j) = (1/(2*pi*var_1_hat))^(T/2)*exp(-T/2);
    end
    t_0_hat = find(L_val == max(L_val));
    t_0_hat_val(i) = t_0_hat;
    mu_0_hat = mu_0_hat_val(t_0_hat);
    delta_hat = delta_hat_val(t_0_hat);
    var_1_hat = var_1_hat_val(t_0_hat);
    % test statistic
    lambda_x = (var_1_hat/var_0_hat)^(T/2);
    chi_statistic = -2*log(lambda_x);
    S(i) = chi_statistic;
end

%% Now try delta = 0 1000 times to investigate the false positive rate. 

delta = 0;
S = zeros(1, N); % empty vector to store test statistics 
t_0_hat_val = zeros(1, N); % empty vector to store t_0_hat values 

for i = 1:N
    X = zeros(1, T);
    epsilon = normrnd(mu_epsilon, sigma, 1, T);
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    var_0_hat = var(X, 1);
    % create a vector of t_0 values, calculate mu_0_hat, delta_hat, var_hat_1
    t_0_val = 1:200; % vector of t_0 values to test
    % initialize vectors to store values of mu_0_hat, delta_hat, var_hat
    % and L (likelihood) that correspond to each t_0 value
    mu_0_hat_val = zeros(1, length(t_0_val));
    delta_hat_val = zeros(1, length(t_0_val));
    var_1_hat_val = zeros(1, length(t_0_val));
    L_val = zeros(1, length(t_0_val));
    for j = 1:length(t_0_val)
        mu_0_hat = mean(X(1:(t_0_val(j) - 1)));
        mu_0_hat_val(j) = mu_0_hat;
        xbar_1_hat = mean(X(1:(t_0_val(j) - 1)));
        xbar_2_hat = mean(X(t_0_val(j):T));
        delta_hat = xbar_2_hat - xbar_1_hat;
        delta_hat_val(j) = delta_hat;
        var_1_hat = (1/T)*((t_0_val(j) - 1)*var(X(1:(t_0_val(j) - 1)), 1) + (T - t_0_val(j) + 1)*var(X(t_0_val(j):T), 1));
        var_1_hat_val(j) = var_1_hat;
        L_val(j) = (1/(2*pi*var_1_hat))^(T/2)*exp(-T/2);
    end
    t_0_hat = find(L_val == max(L_val));
    t_0_hat_val(i) = t_0_hat;
    mu_0_hat = mu_0_hat_val(t_0_hat);
    delta_hat = delta_hat_val(t_0_hat);
    var_1_hat = var_1_hat_val(t_0_hat);
    % test statistic
    lambda_x = (var_1_hat/var_0_hat)^(T/2);
    chi_statistic = -2*log(lambda_x);
    S(i) = chi_statistic;
end

% Calculate the false positive rate of the test 
false_pos_rate = (sum(S > critical_val))/length(S);

% False positive rate is greater than 0.05 = 0.5190

%% Histogram of t_0_hat values
figure()
hist = histogram(t_0_hat_val, 'BinWidth', 1);
title('Estimated t_0 Values', 'FontSize', 14)
xlabel('t_0', 'FontSize', 16)
ylabel('Frequency', 'FontSize', 16)
set(gca,'FontSize', 15) % change font size of axis numbers
saveas(hist, 't_0_histogram', 'png')

%% Delta = 10 

delta = 10;
N = 1000;
S = zeros(1, N);

for i = 1:1000
    X = zeros(1, T);
    epsilon = normrnd(mu_epsilon, sigma, 1, T);
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    var_hat_0 = var(X, 1);
    % create a vector of t_0 values, calculate mu_0_hat, delta_hat, var_hat_1
    t_0_val = 2:200; % vector of t_0 values to test
    % initialize vectors to store values of mu_0_hat, delta_hat, var_hat
    % and L (likelihood) that correspond to each t_0 value
    mu_0_hat_val = zeros(1, length(t_0_val));
    delta_hat_val = zeros(1, length(t_0_val));
    var_1_hat_val = zeros(1, length(t_0_val));
    L_val = zeros(1, length(t_0_val));
    for j = 1:length(t_0_val)
        mu_0_hat = mean(X(1:t_0_val(j) - 1));
        mu_0_hat_val(j) = mu_0_hat;
        delta_hat = mean(X(t_0_val(j):T)) - mean(X(1:t_0_val(j) - 1));
        delta_hat_val(j) = delta_hat;
        var_hat_1 = (1/T)*((t_0_val(j) - 1)*var(X(1:t_0_val(j) - 1), 1) + (T - t_0_val(j) + 1)*var(X(t_0_val(j):T), 1));
        var_1_hat_val(j) = var_hat_1;
        L_val(j) = (1/(2*pi*var_hat_1))^(T/2)*exp((-T/2));
    end
    t_0_hat = find(L_val == max(L_val));
    mu_0_hat = mu_0_hat_val(t_0_hat);
    delta_hat_l = delta_hat_val(t_0_hat);
    var_hat_1_l = var_1_hat_val(t_0_hat);
    % test statistic
    lambda_x = (var_hat_1_l/var_hat_0)^(T/2);
    chi_statistic = -2*log(lambda_x);
    S(i) = chi_statistic;
end

%% Calculate power of the test 
false_neg_rate = sum(S < critical_val)/length(S);
power = 1 - false_neg_rate;

%% Delta vs power
delta = 0.01:0.01:5;  % set delta vector
N = 1000;
obs_delta_t0unknown = zeros(length(delta), N); % empty matrix to store data

for i = 1:length(delta)
    S = zeros(1, N); % empty vector to store chi statistics 
    for j = 1:N
        X = zeros(1, T); % empty vector to store observed data 
        epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 
        for t = 1:T
            X(t) = mu_0 + delta(i)*heaviside(t - t_0) + epsilon(t);
        end
        % create a vector of t_0 values, calculate mu_0_hat, delta_hat, var_hat_1
        t_0_val = 1:200; % vector of t_0 values to test
        % initialize vectors to store values of mu_0_hat, delta_hat, var_hat
        % and L (likelihood) that correspond to each t_0 value
        mu_0_hat_val = zeros(1, length(t_0_val));
        delta_hat_val = zeros(1, length(t_0_val));
        var_1_hat_val = zeros(1, length(t_0_val));
        L_val = zeros(1, length(t_0_val));
        for k = 1:length(t_0_val)
            mu_0_hat = mean(X(1:(t_0_val(k) - 1)));
            mu_0_hat_val(k) = mu_0_hat;
            xbar_1_hat = mean(X(1:(t_0_val(k) - 1)));
            xbar_2_hat = mean(X(t_0_val(k):T));
            delta_hat = xbar_2_hat - xbar_1_hat;
            delta_hat_val(k) = delta_hat;
            var_1_hat = (1/T)*((t_0_val(k) - 1)*var(X(1:(t_0_val(k) - 1)), 1) + (T - t_0_val(k) + 1)*var(X(t_0_val(k):T), 1));
            var_1_hat_val(k) = var_1_hat;
            L_val(k) = (1/(2*pi*var_1_hat))^(T/2)*exp(-T/2);
        end
        t_0_hat = find(L_val == max(L_val));
        if length(t_0_hat) > 1
            t_0_hat = t_0_hat(1);
        end
        mu_0_hat = mu_0_hat_val(t_0_hat);
        delta_hat = delta_hat_val(t_0_hat);
        var_1_hat = var_1_hat_val(t_0_hat);
        var_0_hat = var(X, 1);
        % test statistic
        lambda_x = (var_1_hat/var_0_hat)^(T/2);
        chi_statistic = -2*log(lambda_x);
        S(j) = chi_statistic;
    end
        obs_delta_t0unknown(i, :) = S;
end

%% Plot power of the hypothesis test as function of delta
% (Recall power of a test: probability that test rejects H0 when H1 is true)
% power = 1 - false negative rate 

false_neg_rate_t0_unknown = sum(obs_delta_t0unknown < critical_val, 2)/N;
power = 1 - false_neg_rate_t0_unknown;

figure()
plot = plot(delta(1:100), power(1:100), 'b', 'LineWidth', 1.5); % plot from delta = 0 to 1 
title('t_0 Unknown: Delta vs Power', 'FontSize', 14)
xlabel('Delta', 'FontSize', 16)
ylabel('Power', 'FontSize', 16)
set(gca,'FontSize', 15) % change font size of axis numbers
saveas(plot, 't0_unknown_delta_vs_power', 'png')

