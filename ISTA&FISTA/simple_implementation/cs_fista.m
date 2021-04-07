function [x_2,error] = cs_fista(y,A,lambda,epsilon,itermax)
% Fast Iterative Soft Thresholding Algorithm(FISTA)
% Version: 1.0 written by yfzhang @2019-12-8
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

N = size(A,2);
error = [] ;

x_0 = zeros(N,1);
x_1 = zeros(N,1);
t_0 = 1 ;

for i = 1:itermax
    t_1 = (1+sqrt(1+4*t_0^2))/2 ;
    % g_1 = A'*(y-A*x_1);
    alpha =1;
    % alpha = (g_1'*g_1)/((A*g_1)'*(A*g_1)) ;
    z_2 = x_1 + ((t_0-1)/(t_1))*(x_1 - x_0) ;
    z_2 = z_2+A'*(y-A*z_2);
    x_2 = sign(z_2).*max(abs(z_2)-alpha*lambda,0) ;
    error(i,1) = norm(x_2 - x_1)/norm(x_2) ;
    error(i,2) = norm(y-A*x_2) ;
    if error(i,1) < epsilon || error(i,2) < epsilon
        break;
    else
        x_0 = x_1 ;
        x_1 = x_2 ;
        t_0 = t_1 ;
    end
end
