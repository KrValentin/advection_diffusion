function F=Advdiff_bc_rhs(coeff,lbc,rbc,h,F)

% ADVDIFF_BC_RHS - Add  Dirichlet and/or Robin b.c. in the Right-hand side F
%       for advection-diffusion equation with UPWIND HYBRID scheme 
%------------------------------
% Creation : Caroline Japhet
% Last modification : 31/10/23
%------------------------------
   
b=coeff.b; nu=coeff.nu; 
bp=max(b,0); bm=min(b,0); 
d=2*nu/h;

% ---- left border
if ( strcmp(lbc.type,'Dirichlet') )       
  F(1)=F(1)+(d+bp)*lbc.g;
elseif ( strcmp(lbc.type,'Robin') ) 
  F(1)=F(1)+(d+bp)/(d+bp+lbc.alpha)*lbc.g;
end

% ---- right border
if ( strcmp(rbc.type,'Dirichlet') )
  F(end)=F(end)+(d-bm)*rbc.g;
elseif ( strcmp(rbc.type,'Robin') )
  F(end)=F(end)+(d-bm)/(d-bm+rbc.alpha)*rbc.g;
end
 
