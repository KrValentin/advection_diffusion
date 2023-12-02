function [X,h,t,dt,lg,rg,U0]=Meshes_data(xa,xb,nx,t0,tf,nt,lbc,rbc,u0);
%
% MESHES_DATA - space-time mesh and datas (b.c. and u0) on the meshes
%-----------------------------------------------------------------------------
% Creation : Caroline Japhet
% Last modification : 31/10/23
%-----------------------------------------------------------------------------

T =linspace(xa,xb,nx+1);                                     % meshes T and D
h=T(2)-T(1);
X=[xa T(1:end-1)+h/2 xb]';
t=linspace(t0,tf,nt+1)';  dt=(tf-t0)/nt;                     % time grid  

lg=lbc.g(xa,t(2:end)); rg=rbc.g(xb,t(2:end));                % b.c. on gam
U0=u0(X);                                                    % initial condition
