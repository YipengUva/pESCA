%% test fast SVD

%generate data sets
n = 1000;
p = 300;
X = randn(n,p);
Y = randn(p,n);

% standard trace function
tic
for i = 1:10,
    [result] = trace(X*Y);
end
trace_time = toc;

% fast trace
tic
for i = 1:10,
    [result_fast] = fastTrace(X,Y);
end
fastTrace_time = toc;

% test the results
result
result_fast

% compare the time 
trace_time
fastTrace_time


