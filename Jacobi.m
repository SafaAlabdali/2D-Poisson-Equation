%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using diagonal matrix DT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
%create matrix
N=4; 
alpha=0.05;
beta =-0.6;
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
D=diag(diag(A)); % diagonal matrix 
u0=zeros(n,1); % initial vector 

k=180; % number of iteration 
r=zeros(k,1);  % residual 
e=zeros(k,1);  % error 

for i=1:k
    %uk=u0 + inv(D) *(b-A*u0)
    uk=u0+( D \(b-A*u0) );
    r(i)=norm(A*uk-b);
    e(i)=norm(uk-u) ;
    u0=uk;
end 

spy(A)
title('matlab spy plot for matrix A')
figure 
semilogy(1:k,r);
hold on 
semilogy(1:k,e);
legend ('Residual','Error');
title('Jacobi Method with matrix D')