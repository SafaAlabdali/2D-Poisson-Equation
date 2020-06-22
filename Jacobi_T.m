%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using tridiagonal matrix T 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
%create matrix
N=10; 
alpha=1;
beta =1;
dx=0.5;
dy=0.4;
A=A2D(N,alpha,beta,dx,dy);

n=size(A);
n=n(1);
% create random vector u 
u=rand(n,1);
% find the vecto b 
b=A*u;
%solve for Au=b
% construct tridiagonal mat 
T = zeros(n,n);
T(1:1+n:n*n) = diag(A);
T(n+1:1+n:n*n) = diag(A,1);
%T(2:1+N:N*N-N) = diag(A,-1);
u0=zeros(n,1); % initial vector 
r=zeros(n,1);  % residual 
e=zeros(n,1);  % error 

k=1000; % number of iteration 
for i=1:k 
    r=b-A*u0;
    uk=u0+ solve(T,r);
    r(i)=norm(A*uk-b);
    e(i)=norm(uk-u) ;
    u0=uk;
end 

spy(A)
title('matlab spy plot for matrix A')
figure 
loglog(1:k,r);
hold on 
loglog(1:k,e);
legend ('Residual','Error');
title('Jacobi Method with matrix T')













