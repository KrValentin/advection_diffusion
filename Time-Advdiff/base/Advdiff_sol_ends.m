function U=Advdiff_sol_ends(coeff,lbc,rbc,h,U)

% ADVDIFF_SOL_ENDS - adds in U the values at ends of the domain
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
  U0=lbc.g;
elseif ( strcmp(lbc.type,'Robin') )      
  U0=(lbc.g +(d-bm)*U(1))/(d+bp+lbc.alpha);
end
  
% ---- right border
if ( strcmp(rbc.type,'Dirichlet') )
  Unp1=rbc.g;
elseif ( strcmp(rbc.type,'Robin') )      
  Unp1=(rbc.g +(d+bp)*U(end))/(d-bm+rbc.alpha);
end

U=[U0; U; Unp1];  
