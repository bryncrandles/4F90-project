% Let us just try to create the permutation test using delta = 5 to begin

% Generate time series data 
T = 200; % length of time 
t_0 = 100; % time of change point 
mu_0 = 10;
N = 1000; % number of simulations 
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error
delta = 5;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end

%% Calculate the cusum statistic 
% first calculate mean
xbar = mean(X);
% immediately calculate the normalized cusum statistics S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, (T - 1));
for t = 1:(T - 1)
    S(t) = sum(X(1:t) - xbar);
end
% normalize
S_N = zeros(1, (T - 1));
for t = 1:(T - 1)
    S_N(t) = sqrt(T/(t*(T - t)))*abs(S(t));
end

%% immediately calcualte normalized - check that it is the same

S_test = zeros(1, (T - 1));
for t = 1:(T - 1)
    S_test(t) = sqrt(T/(t*(T - t)))*abs(sum(X(1:t) - xbar));
end
% yes they are the same 

%% plot

figure()
plot(1:(T - 1), S_test);

%% the cusum statistic is the maximum of all of these 
c_stat = max(S_N);
% estimated value of t_0
t_0_hat = find(S_N == c_stat);

%% permute X(t) and calculate cusum statistic for that permutation
X_rand = X(randperm(T));

% calculate statistic
S_N_rand = zeros(1, (T - 1));
% loop to calculate S(t)
for t = 1:(T - 1)
    S_N_rand(t) = sqrt(T/(t*(T - t)))*abs(sum(X_rand(1:t) - xbar));
end

c_stat_rand = max(S_N_rand);

%% do this 100000 times
P = 10000;
c_stat_rand = zeros(1, P);
for i = 1:P
    X_rand = X(randperm(T));
    S_N_rand = zeros(1, (T - 1));
    % loop to calculate S(t)
    for t = 1:(T - 1)
        S_N_rand(t) = sqrt(T/(t*(T - t)))*abs(sum(X_rand(1:t) - xbar));
    end
    c_stat_rand(i) = max(S_N_rand);
end

%%
figure()
histogram(c_stat_rand, 'BinWidth', 0.1)
p_val = sum(c_stat_rand > c_stat)/length(c_stat_rand);


