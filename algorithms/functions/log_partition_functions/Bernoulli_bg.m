function [ out ] = Bernoulli_bg(Theta)
% The first order gradient of the log-partation function of Bernoulli distribution;
% The following is the matrix form
out = exp(Theta)./(1+exp(Theta));
end

