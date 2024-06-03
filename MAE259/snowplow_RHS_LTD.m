function f=snowplow_RHS_LTD(x,t,p)

dt=1e-11;
%[L,q,I,M,r,v,Is]
L  = real(x(1)); % load inductance
Q  = real(x(2)); % charge on capacitor
I  = real(x(3)); % load current
M  = real(x(4)); % piston mass
R  = real(x(5)); % piston radius
V  = real(x(6)); % piston velocity
Is = real(x(7)); % shunt current

rho_gain = 2*pi*R*p.mu*p.mp*(snowplow_density_LTD(R,p)+snowplow_density_LTD(R+V*dt,p))/2.0;

f=zeros(size(x));

f(1)=-p.mu0*p.l/(2*pi)*V/R;

f(2)=I+Is;

f(3)=(Is*p.Rsh - I*(p.Rpf+f(1)))/(p.Lpf+L);

f(4)= -rho_gain*V;

f(5)=V;

f(6)=(-f(4)*V-I+p.P_coeff*(p.Rf/10/R)^4*((p.Rf/R)^(7/3)-1))/M;

f(7)= (-Q/p.Cbr - Is*p.Rsh - (I+Is)*p.Rbr - p.Lbr*f(3))/p.Lbr;

end
