% EX_ADVDIFF : Solve time-dependent advection-diffusion equation with Dirichlet and/or Robin b.c. :
%
% du/dt + eta*u  + b*du/dx - nu*d^2u/dx^2 = f   in [X(1),X(end)]x(t0,tf)
%
% IF lbc.type = 'Dirichlet'
% u(X(1),.) = lbc.g(.) on (t0,tf)
% IF lbc.type = 'Robin'
% (-nu*u'+b+lbc.alpha*u)(X(1),.)=lbc.g(.) on (t0,tf)
%
% IF rbc.type = 'Dirichlet'
% u(X(end),.) = rbc.g(.) on (t0,tf)
% IF rbc.type = 'Robin'
% (nu*u'-b+rbc.alpha*u)(X(end),.)=rbc.g(.) on (t0,tf)
%
% with a "cell-centered" finite volume method on a uniform mesh			  	
% X is an (n+2)x1 array, representing mesh D, of [X(1),X(end)]
% (with mesh size h/2 on first and last intervals, and mesh size h else, with h=1/n).
%
%   |-----------|-----------| ...  |-----------|      mesh T
%   o-----o-----|-----o-----| ...  |-----o-----o      mesh D
%  x0     x1          x2                 xn    xn+1
%
% UPWIND HYBRID scheme for advection, of order 1 (and order 2 for diffusion)
% and implicite Euler scheme in time (order 1)
%------------------------------
% Creation : Caroline Japhet
% Last modification : 31/10/23
%------------------------------

addpath data:base:graphics
format short e

%[coeff,f,lbc,rbc,xa,xb,nx,t0,tf,nt,u0,uex,PLOTSOL]=Data_test_scheme();  % data
[coeff,f,lbc,rbc,xa,xb,nx,t0,tf,nt,u0,uex,PLOTSOL]=Data();  

[X,h,t,dt,lbc.g,rbc.g,U0]=Meshes_data(xa,xb,nx,t0,tf,nt,lbc,rbc,u0);

Uex=zeros(length(X),nt+1);                                 
for n=1:nt+1
  Uex(:,n)=uex(X,t(n));                                     % exact solution
end

A=Advdiff_matrix(coeff,lbc,rbc,nx,h,dt);
F=Advdiff_rhs(f,X,h,t);
U=Advdiff(coeff,h,nt,dt,lbc,rbc,A,F,U0);                    % discrete solution

n=nt+1;
Plotsol(X,U(:,n),Uex(:,n),PLOTSOL);                         % plots

[errLinfL2,errLinf]=Errsol(U,Uex);                              % L2 and Linf errors 

