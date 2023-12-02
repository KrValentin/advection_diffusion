function [coeff,f,lbc,rbc,xa,xb,nx,t0,tf,nt,u0,u,PLOTSOL]=Data_test_scheme();

% DATA_TEST_SCHEME - datas for testing the VF scheme :
% advection : the scheme is exact with a polynomial of degree 0 in space
% diffusion : the scheme is exact with a polynomial of degree <= 1 in space
%-----------------------------------------------------------------------------
% Creation : Caroline Japhet
% Last modification : 24/10/23
%-----------------------------------------------------------------------------

% --- CHOOSE FOR THE TEST ---------------------  
DEG=0    % '0' for u=2*t+3 and '1' for u=(3x+1)(2*t+3)
test=1   % '0' for Dirichlet and '1' for Robin
%----------------------------------------------

% ---- physical data --------------------------
xa = 0; xb = 1;   % global space domain [xa,xb] 
t0 = 0; tf = 1;   % global time domain [t0,tf] 
eta=1;            % reaction 
b=1;              % advection
nu=0.1;           % diffusion
coeff.eta=eta; coeff.b=b; coeff.nu=nu;
%----------------------------------------------

% ---- Numerical data -------------------------
nx=5;     % number of intervals of mesh T
nt=5;    % number of time steps
%----------------------------------------------

% ---- exact solution -------------------------
if DEG==0
  u=@(x,t) ones(size(x)).*(2*t+3);
  ut=@(x,t) 2*ones(length(x),length(t));
  ux=@(x,t) zeros(length(x),length(t));
  uxx=@(x,t) zeros(length(x),length(t));
elseif DEG==1
  u=@(x,t) (3*x+1).*(2*t+3);
  ut=@(x,t) 2*(3*x+1).*ones(size(t));
  ux=@(x,t) 3*(2*t+3);
  uxx=@(x,t) zeros(length(x),length(t));
end
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



