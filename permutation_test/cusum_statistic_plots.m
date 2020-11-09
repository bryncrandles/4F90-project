% Generate time series data 

T = 200; % length of time 
t_0 = 100; % time of change point 
mu_0 = 10;
N = 1000; % number of simulations 
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error

%% delta = 0
delta = 0;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end
% calculate mean
xbar = mean(X);
% cumulative summation time series S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, T);
for t = 1:T
    S(t) = sum(X(1:t) - xbar);
end

figure()
plot(1:T, S)

%%
% normalize
S_N = zeros(1, T);

for t = 1:T
    S_N(t) = sqrt(T/(t*(T - t)))*abs(S(t));
end

figure()
plot(1:T, S_N)

%% try with delta = 10 and t0 = 100
delta = 10;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end
% calculate mean
xbar = mean(X);
% cumulative summation time series S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, T);
for t = 1:T
    S(t) = sum(X(1:t) - xbar);
end

figure()
plot(1:T, S)

% normalize
S_N = zeros(1, T);

for t = 1:T
    S_N(t) = sqrt((T)/(t*(T - t)))*abs(S(t));
end

figure()
plot(1:T, S_N)

%% try t0 = 50 with delta = 0
t_0 = 50;
delta = 0;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end
% calculate mean
xbar = mean(X);
% cumulative summation time series S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, T);
for t = 1:T
    S(t) = sum(X(1:t) - xbar);
end

figure()
plot(1:T, S)

% normalize
S_N = zeros(1, T);

for t = 1:T
    S_N(t) = sqrt(T)/(t*(T - t))*abs(S(t));
end

figure()
plot(1:T, S_N)

%% try t0 = 50 with delta = 10
t_0 = 50;
delta = 10;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end
% calculate mean
xbar = mean(X);
% cumulative summation time series S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, T);
for t = 1:T
    S(t) = sum(X(1:t) - xbar);
end

figure()
plot(1:T, S)

% normalize
S_N = zeros(1, T);

for t = 1:T
    S_N(t) = sqrt(T/(t*(T - t)))*abs(S(t));
end

figure()
plot(1:T, S_N)

%% try t0 = 150 with delta = 0
t_0 = 150;
delta = 0;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end
% calculate mean
xbar = mean(X);
% cumulative summation time series S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, T);
for t = 1:T
    S(t) = sum(X(1:t) - xbar);
end

figure()
plot(1:T, S)

% normalize
S_N = zeros(1, T);

for t = 1:T
    S_N(t) = sqrt(T)/(t*(T - t))*abs(S(t));
end

figure()
plot(1:T, S_N)
t_0 = 50;
delta = 0;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end
% calculate mean
xbar = mean(X);
% cumulative summation time series S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, T);
for t = 1:T
    S(t) = sum(X(1:t) - xbar);
end

figure()
plot(1:T, S)

% normalize
S_N = zeros(1, T);

for t = 1:T
    S_N(t) = sqrt(T)/(t*(T - t))*abs(S(t));
end

figure()
plot(1:T, S_N)

%% try t0 = 150 with delta = 10

t_0 = 150;
delta = 10;
X = zeros(1, T); % empty vector to store observed data 
epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 

% loop to generate time series 
for t = 1:T
     X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
end
% calculate mean
xbar = mean(X);
% cumulative summation time series S(t) = sum(Xi, i = 1..t) - xbar
S = zeros(1, T);
for t = 1:T
    S(t) = sum(X(1:t) - xbar);
end

figure()
plot(1:T, S)

% normalize
S_N = zeros(1, T);

for t = 1:T
    S_N(t) = sqrt(T)/(t*(T - t))*abs(S(t));
end

figure()
plot(1:T, S_N)




    