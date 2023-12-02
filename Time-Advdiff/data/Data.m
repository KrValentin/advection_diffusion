function [coeff,f,lbc,rbc,xa,xb,nx,t0,tf,nt,u0,u,PLOTSOL]=Data();

% DATA - datas for testing the implicit Euler / VF scheme and for the DD code
% advection : the scheme is of order 1 in space
% diffusion : the scheme is of order 2 in space
% time scheme : order 1
%-----------------------------------------------------------------------------
% Creation : Caroline Japhet
% Last modification : 24/10/23
%-----------------------------------------------------------------------------

% --- CHOOSE FOR THE TEST ---------------------  
test=1   % '0' for Dirichlet and '1' for Robin
%----------------------------------------------

% ---- physical data --------------------------
xa = 0; xb = 1;   % global space domain [xa,xb] 
t0 = 0; tf = 1;   % global time domain [t0,tf] 
eta=0;            % reaction 
b=1;              % advection
nu=0.1;           % diffusion
coeff.eta=eta; coeff.b=b; coeff.nu=nu;
%----------------------------------------------

% ---- Numerical data -------------------------
nx=100;     % number of intervals of mesh T
nt=100;    % number of time steps
%----------------------------------------------

% ---- exact solution -------------------------
u=@(x,t) sin(3*x).*cos(2*t);
ut=@(x,t) -2*sin(3*x).*sin(2*t);
ux=@(x,t) 3*cos(3*x).*cos(2*t);
uxx=@(x,t) -9*sin(3*x).*cos(2*t);
%----------------------------------------------

u0 = @(x) u(x,t0);
f = @(x,t)  ut(x,t)+eta*u(x,t)+b*ux(x,t)-nu*uxx(x,t);  % source term

if test==0
  lbc.type='Dirichlet'; 
  rbc.type='Dirichlet';
  lbc.g=@(x,t) u(x,t);
  rbc.g=@(x,t) u(x,t);
elseif test==1
  lbc.type='Robin';
  rbc.type='Robin';
  pa=1; pb=2;
  lbc.alpha=pa-b/2;         % Robin parameters
  rbc.alpha=pb+b/2;       
  lbc.g=@(x,t) -nu*ux(x,t)+(b+lbc.alpha)*u(x,t);
  rbc.g=@(x,t) nu*ux(x,t)+(-b+rbc.alpha)*u(x,t);
end



% ---- For graphics ----------------------
PLOTSOL=1;  % '1' plot solutions, '0' else
%-----------------------------------------   



