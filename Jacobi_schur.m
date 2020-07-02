%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using tridiagonal matrix T 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
%create matrix
N=16; 
alpha=1;
beta =1E4;
Lx=1.0; % Physical size of the domain in X-direction
Ly=0.1; % Physical size of the domain in Y-direction
dx=Lx/N;
dy=Ly/N;
A=A2D(N,alpha,beta,dx,dy);

n=size(A);
n=n(1);
% create random vector u 
u=rand(n,1);
% find the vecto b 
b=A*u;
%solve for Au=b
% construct tridiagonal mat 
D = zeros(n,n);
D(1:1+n:n*n) = diag(A);
T = zeros(n,n);
T(1:1+n:n*n) = diag(A);
T(n+1:1+n:n*n) = diag(A,1);
%T(2:1+N:N*N-N) = diag(A,-1);
u0=zeros(n,1); % initial vector 
k=100; % number of iteration 
r=zeros(n,1);  % residual 
e=zeros(1,k);  % norm of error at every iteration
res=zeros(1,k); % norm of residual at every iteration
res_0=norm(A*u0-b);
e_0=norm(u0-u) ;
for i=1:k 
    r=b-A*u0;
    uk=u0+ solve(A,r);
    res(i)=norm(A*uk-b);
    e(i)=norm(uk-u) ;
    u0=uk;
end 

spy(A)
title('matlab spy plot for matrix A')
figure 
semilogy(1:k,res/res_0);
hold on 
semilogy(1:k,e/e_0);
legend ('Residual','Error');
title('Jacobi Method with Schur-complement approach')