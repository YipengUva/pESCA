function [ out ] = Bernoulli_b(Theta)
% The log-partation function of Bernoulli distribution;
% The following is the matrix form.

out = log(1+exp(Theta));

end

