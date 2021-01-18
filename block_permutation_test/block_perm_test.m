% block permutation test for correlated errors using AR(3) process

B = 5; % number of blocks
T = 200; % length of time series
W = T/B; % number of windows

t_0 = 100; % time of change point 
delta = 10;
N = 1000; % number of simulations 
phi0 = 0;%5.5; % set to 5.5 for mean of 10
phi1 = 0.2; 
phi2 = 0.15;
phi3 = 0.1;
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error
X = zeros(1, T);
X(1) = 0; 
X(2) = 0;
X(3) = 0;
P = 10000;
stat = zeros(1, N);
alpha = 0.05;

for i = 1:N
    epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error
    for t = 4:T
        X(t) = phi0 + delta*heaviside(t - t_0) + phi1*X(t - 1) + phi2*X(t - 2) + phi3*X(t - 3) + epsilon(t);
    end    
    xbar = mean(X);
    S = zeros(1, (T - 1));
    for t = 1:(T - 1)
        S(t) = sqrt(T/(t*(T - t)))*abs(sum(X(1:t) - xbar));
    end
    c_stat = max(S);
    c_stat_rand = zeros(1, P);
    c_stat_rand(1) = c_stat;
    for j = 2:P
        R = reshape(X, W, B);
        R_rand = R(:, randperm(B));
        X_rand = reshape(R_rand, 1, T);
        S_N_rand = zeros(1, (T - 1));
        % loop to calculate S(t)
        for t = 1:(T - 1)
            S_N_rand(t) = sqrt(T/(t*(T - t)))*abs(sum(X_rand(1:t) - xbar));
        end
        c_stat_rand(j) = max(S_N_rand);
    end
    p_val = sum(c_stat_rand >= c_stat)/P;
    stat(i) = p_val;
end

false_neg_rate = sum(stat > alpha, 2)/N;
power = 1 - false_neg_rate;

save('block_perm_test_power.mat', 'power', '-mat')

exit

