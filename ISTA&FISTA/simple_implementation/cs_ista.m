function [x_hat,error] = cs_ista(y,A,lambda,epsilon,itermax)
% Iterative Soft Thresholding Algorithm(ISTA)
% Version: 1.0 written by Louis Zhang @2019-12-7
% Reference: Beck, Amir, and Marc Teboulle. "A fast iterative 
% shrinkage-thresholding algorithm for linear inverse problems." 
% SIAM journal on imaging sciences 2.1 (2009): 183-202.

% Inputs:
% y         - measurement vector
% A         - measurement matrix
% lambda    - denoiser parameter in the noisy case
% epsilon   - error threshold
% inter_max - maximum number of amp iterations
%
% Outputs:
% x_hat     - the last estimate
% error     - reconstruction error

if nargin < 5
    itermax = 10000 ;
end
if nargin < 4
    epsilon = 1e-4 ;
end
if nargin < 3
    lambda = 2e-5 ;
end

N = size(A,2) ;
error = [];
x_1 = zeros(N,1) ;

for i = 1:itermax
    g_1 = A'*(y - A*x_1) ;
    alpha = 1 ;
    % obtain step size alpha by line search
    % alpha = (g_1'*g_1)/((A*g_1)'*(A*g_1)) ;
    x_2 = x_1 + alpha * g_1 ;
    x_hat = sign(x_2).*max(abs(x_2)-alpha*lambda,0) ;
    error(i,1) = norm(x_hat - x_1) / norm(x_hat) ;
    error(i,2) = norm(y-A*x_hat) ;
    if error(i,1) < epsilon || error(i,1) < epsilon
        break;
    else
        x_1 = x_hat ;
    end
end
