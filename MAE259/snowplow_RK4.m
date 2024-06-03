function [xf,tf]=RK4(R,x,t,s,p,plotter)

h=s.h; xold=x; ooo=0;

for n=1:s.T/h;
  if p.config==1, radius=x(5);, initial_radius=p.Rf;, end
  if p.config==0,radius=x(2);,initial_radius=p.r0;, end
  if initial_radius/radius>s.CR, break, end
  f1=feval(R,x,t,p); f2=feval(R,x+h*f1/2,t,p); f3=feval(R,x+h*f2/2,t,p); f4=feval(R,x+h*f3,t,p);
  x=x+h*((f1+f4)/6+(f2+f3)/3);
  t=t+h;
  xxx=feval(@snowplow_Vg,t);
  %feval(plotter,xold,x,t-h,t,ooo,xxx); xold=x; ooo=xxx;
  xf(:,n)=x;
  tf(n)=t;
end
end
