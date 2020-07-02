%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using diagonal matrix DT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % to make it clear 
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
D=diag(diag(A)); % diagonal matrix 
u0=zeros(n,1); % initial vector 
r=zeros(n,1); % residual

k=100; % number of iterations
res=zeros(k,1);  % residual 
e=zeros(k,1);  % error 
res_0=norm(A*u0-b);
e_0=norm(u0-u) ;
for i=1:k
    %uk=u0 + inv(D) *(b-A*u0)
    r=b-A*u0;
    uk=u0+ D\r;
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
title('Jacobi Method with matrix D')