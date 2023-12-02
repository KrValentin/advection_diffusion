function A=Advdiff_matrix(coeff,lbc,rbc,nx,h,dt)

% ADVDIFF_MATRIX - computes matrix with Dirichlet and/or Robin b.c.
%       for advection-diffusion equation with UPWIND HYBRID scheme 
%------------------------------
% Creation : Caroline Japhet
% Last modification : 31/10/23
%------------------------------
   
eta=coeff.eta; b=coeff.b; nu=coeff.nu;
eta=eta+1/dt;
bp=max(b,0); bm=min(b,0); ba=abs(b);
d=2*nu/h;

% global matrix (with sparse format)
e=ones(nx,1);  
Cm=(d-bm)^2/(2*d+ba); Cp=(d+bp)^2/(4*nu/h+ba); 
  
A=spdiags([-Cp*e,(eta*h+Cm+Cp)*e,-Cm*e],[-1,0,1],nx,nx);

% add b.c. in A 			  
%-----------------------
% left border
if ( strcmp(lbc.type,'Dirichlet') )       
  A(1,1)=eta*h+Cm+(d+bp);
elseif ( strcmp(lbc.type,'Robin') ) 
  A(1,1)=eta*h+Cm+(d+bp)*(lbc.alpha+b)/(d+bp+lbc.alpha);
end

% right border
if ( strcmp(rbc.type,'Dirichlet') )
  A(nx,nx)=eta*h+Cp+(d-bm);
elseif ( strcmp(rbc.type,'Robin') )      
  A(nx,nx)=eta*h+Cp+(d-bm)*(rbc.alpha-b)/(d-bm+rbc.alpha); 
end
   
