function dxdt = eqn_RLCcore_Zpinch_snowplow(t,x,params,dens_prof,dens_par) %#ok<INUSL>
%#ok<*NASGU>
Lbr = params(1);  %Brick inductance
Rbr = params(2);  %Brick resistance
Cbr = params(3);  %Circuit capacitance
Lpf = params(4);  %Power feed inductance
Rpf = params(5);  %Power feed resistance
Rsh = params(6);  %Shunt resistance
mu = params(7);   %Atomic mass
mml = params(8);  %Total sheath mass
gam = params(9);  %Gamma
Rext = params(10);     %Return current radius
Hext = params(11);     %Pinch height
dt = params(12);       %Timestep
P_coeff = params(13); % Pressure coefficient
R_max = params(14);   %Initial radius
R_min = params(15);   %Minimum radius

%% physical constants
c = 299792458;
mp = 1.0e-3/6.02214086e23;
mu0 = 4*pi*1.0e-7;
Ev_J = 1.602e-19;

%% Define the variables to make things little more legible
L  = real(x(1)); % load inductance
Q  = real(x(2)); % charge on capacitor
I  = real(x(3)); % load current
M  = real(x(4)); % piston mass
R  = real(x(5)); % piston radius
V  = real(x(6)); % piston velocity
Is = real(x(7)); % shunt current

%% calculate the swept density at this position
rho_gain = 2*pi*R*mu*mp*(dens_prof(R,dens_par)+dens_prof(R+V*dt,dens_par))/2.0;

%% the array dxdt is the same length as x
dxdt = zeros(size(x));

%% define equations for time derivatives

% dL/dt
dxdt(1) = -(mu0*Hext/(2*pi))*V/R;

% dq/dt = I
dxdt(2) = I+Is;

% dI/dt circuit eq. (ORIGINAL)
%dxdt(3) = -1/(Lbr+L)*((Rbr+dxdt(1))*I + Q/Cbr);

% dIload/dt circuit eq. w/ shunt:
% Is*Rsh = I*Rpf + (Lpf+L)*dI/dt + I*dL/dt, solve for dI/dt:
% dI/dt = (Is*Rsh- I*(Rpf+dL/dt))  /(Lpf+L)
dxdt(3) = (Is*Rsh - I*(Rpf+dxdt(1)))/(Lpf+L);

% dM/dt
dxdt(4) = -rho_gain*V;
if dxdt(4)<0.0 || M >= mml % only if mass can be gained
  dxdt(4)=0.0;
end

% dr/dt = v
dxdt(5) = V;

% m dV/dt + V dm/dt = F (units in force/length)
P_coeff2 = P_coeff*(R_min/R)^4;
dxdt(6) = (1/M)*(-dxdt(4)*V - (mu0*I^2)/(4*pi*R) + P_coeff2*((R_max/R)^(7/3)-1) );

% dIshunt/dt
% Vcap = (I+Is)*Rbr + Lbr*(dI/dt+dIs/dt) + Is*Rsh, solve for dIs/dt
% dIs/dt = (-Q/Cbr - Is*Rsh - (I+Is)*Rbr - Lbr*dI/dt)/Lbr
dxdt(7) = (-Q/Cbr - Is*Rsh - (I+Is)*Rbr - Lbr*dxdt(3))/Lbr;
