function f=snowplow_RHS(x,t,p)

dt=1e-9;

%[I,rp,up,M]

I    = real(x(1));
rp   = real(x(2));
up   = real(x(3));
M    = real(x(4));

rho_gain = 2*pi*rp*(snowplow_density(rp,p)+snowplow_density(rp+up*dt,p))/2.0*p.l;

f=zeros(size(x));

f(1) = (feval(@snowplow_Vg,t)-(p.Rg-p.mu0*p.l/(2*pi)*up/rp)*I)/(p.Lg+p.mu0*p.l/(2*pi)*log(p.r0/rp));

f(2) = up;

f(3) = (-p.mu0*I^2/(4*pi*rp)-up*f(4))/M;

f(4) = -rho_gain*up;

%{
f=[(feval(@snowplow_Vg,t)-(p.Rg-2*p.l/p.c^2*x(3)/x(2))*x(1))/(p.Lg+2*p.l/p.c^2*log(p.r0/x(2))),...
x(3),...
(-p.l*x(1)^2/(p.c^2*x(2))+2*pi*x(2)*x(3)^2*feval(@snowplow_density,x(2),p)*p.l)/x(4),...
-2*pi*x(2)*x(3)*feval(@snowplow_density,x(2),p)*p.l
];
%}
end
%6e5
%1e7
