function Plotsol(X,U,Uex,PLOTSOL)

if PLOTSOL==1
  % Plot exact and discrete solutions
  figure(1)
  clf;
  plot(X,U,'m*',X,Uex,'b')
  xlabel('x','Fontsize',14)
  legend('discrete solution','exact solution')
end  
