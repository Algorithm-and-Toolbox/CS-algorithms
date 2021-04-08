% One-dimensional random Gaussian signal test script for CS reconstruction
% algorithm
% Version: 1.0 written by yfzhang @2019-12-8
clear
clc
N = 1024 ;
M = 512 ;
K = 10 ;
x = zeros(N,1);
T = 5*randn(K,1);
index_k = randperm(N);
x(index_k(1:K)) = T;

A = randn(M,N);
A=sqrt(1/M)*A;
A = orth(A')';
% sigma = 1e-4 ;
% e = sigma*randn(M,1);
y = A * x ;% + e ;

[x_rec1,error1] = cs_fista(y,A,5e-3,1e-4,5e3) ;
[x_rec2,error2] = cs_ista(y,A,5e-3,1e-4,5e3) ;

figure
plot(error1(:,2),'r-');
hold on
plot(error2(:,2),'b-');
legend('fista','ista')
hold off

figure
plot(1:1024,x,'r', 1:1024,x_rec1,'g*', 1:1024,x_rec2,'b^');
legend('original','fista','ista')
