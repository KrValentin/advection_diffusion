function F=Advdiff_rhs(f,X,h,t)

% ADVDIFF_RHS - computes physical RHS from source term f, without b.c.
%    for time-dependent advection-diffusion equation with UPWIND HYBRID scheme 
%------------------------------
% Creation : Caroline Japhet
% Last modification : 31/10/23
%------------------------------

nt=length(t)-1;
  
F=zeros(length(X)-2,nt-1);

% -- Time loop --
for n=1:nt
  F(:,n)=h*f(X(2:end-1),t(n+1));     
end
  
