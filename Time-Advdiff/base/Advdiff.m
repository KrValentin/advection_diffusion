function U=Advdiff(coeff,h,nt,dt,lbc,rbc,A,F,U0)

% ADVDIFF : Solve time-dependent advection-diffusion equation with Dirichlet and/or Robin b.c. 
%   with UPWIND HYBRID scheme for advection and implicite Euler scheme in time
%------------------------------
% Creation : Caroline Japhet
% Last modification : 31/10/23
%------------------------------
   
U=U0;    % initial condition
lbcn=lbc; rbcn=rbc;

% -- Time loop --
for n=1:nt
  Fn=h/dt*U(2:end-1,n)+F(:,n);
  lbcn.g=lbc.g(n); rbcn.g=rbc.g(n);
  Fn=Advdiff_bc_rhs(coeff,lbcn,rbcn,h,Fn);             % add B.C. in F
  Un=A\Fn;                                             % solve the system         
  Un=Advdiff_sol_ends(coeff,lbcn,rbcn,h,Un);           % adds in U the values at ends of the domain
  U=[U  Un];
end
