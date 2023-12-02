function [errLinfL2,errLinf]=Errsol(U,V);
%
% ERRSOL : relative Linf-L2 and Linf errors between U and V solutions
%
%-----------------------------------------------------------------------------
% Creation : Caroline Japhet
% Last modification : 31/09/23
%-----------------------------------------------------------------------------

for n=1:size(V,2)
  eL2(n)=norm(U(:,n)-V(:,n));
  nL2(n)=norm(V(:,n));
  eLf(n)=norm(U(:,n)-V(:,n),inf);
  nLf(n)=norm(V(:,n),inf);
end

if norm(V,inf)==0
  errLinfL2=norm(eL2,inf);       
  errLinf=norm(eLf,inf);
else
  errLinfL2=norm(eL2,inf)/norm(nL2,inf);       
  errLinf=norm(eLf,inf)/norm(nLf,inf);
end

disp(' ')
disp('   LinfL2-Error     Linf-Error')
disp([errLinfL2  errLinf])
