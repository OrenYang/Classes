clear all
close all

p.config=1;

%% ALL UNITS IN SI (MKS) except where noted

%% physical constants
p.c = 299792458;
p.mp = 1.0e-3/6.02214086e23;
p.mu0 = 4*pi*1.0e-7;

%% circuit parameters
p.Cbr = 20*(58e-9/2.0);
p.Rbr = 15e-3;   % bricks resistance (parallel)
p.Lbr = 9.0e-9;  % bricks inductance (parallel)
p.Lpf = 10.0e-9; % feed inductance
p.Rpf = 6e-3;    % feed resistance
p.Rsh = 1.9;     % shunt (core) resistance
Vch = 110e3;   % charge voltage

p.Rret = 6.0e-2; % return current radius
p.l = 1.4e-2; % AK gap

%% load parameters
p.mu = 40.0; % atomic mass
p.n0 = 0.0e23*(4.0/p.mu); % peak density - CENTER
p.R0 = 0.0e-2; % initial radius - CENTER
p.DR0 = 4.0e-3; % initial width - CENTER
p.n1 = 3.6e22; % peak density - SHELL
p.R1 = 1.27e-2; % initial radius - SHELL
p.DR1 = 3.5e-3; % initial width - SHELL
p.gamma = 5/3; % polytropic index
p.P_coeff = 1.0e4; % pressure multiplier
p.Rf = p.R1 + 4*p.DR1; %initial radius
Lstart = (p.mu0/(2*pi))*p.l*log(p.Rret/p.Rf); %initial load inductance

%% initial density profile
Rgrid = 0:1e-8:p.Rf; % grid to calculate mass
n = snowplow_density_LTD(Rgrid,p); % initial denstity profile
rho = n*p.mu*p.mp; % initial mass profile
ml = trapz(Rgrid,2*pi.*Rgrid.*rho); % total mass/length
w_bd = 0.5e-3; % initial sheath width
R_bd = p.Rf:1.0e-6:p.Rf+w_bd; % grid to calculate initial mass
N_bd = snowplow_density_LTD(R_bd,p); % initial sheath density
M0 = trapz(R_bd,2*pi*p.mu*p.mp.*R_bd.*N_bd); % initial sheath mass
p.mml=M0+ml;

figure(4), plot(Rgrid,rho,'-b','linewidth',2), hold on, title('Density Profile'),xlabel('Radius [m]'),ylabel('density [kg/m-3]')

t0=0;
x0=[Lstart -Vch*p.Cbr 0.0 M0 p.Rf 0.0 0.0] % initial conditions

s.h=1e-11;
s.T=5e-7;
s.CR=10;

[x,t]=feval(@snowplow_RK4,@snowplow_RHS_LTD,x0,t0,s,p,@snowplow_plotter);

figure(1), plot(t(:),x(5,:),'DisplayName','Radius','linewidth',2), ylim([0,0.03]), hold on, title('Radius vs Time'),legend,xlabel('Time [s]'), ylabel('Radius [m]')
figure(2), plot(t(:),x(3,:),'-r','DisplayName','Current','linewidth',2), ylim([0,6e5]), hold on, title('Current vs Time'),legend,xlabel('Time [s]'), ylabel('Current [A]')
figure(3), plot(t(:),x(6,:),'DisplayName','Velocity'), ylim([-4e5,2e5]), hold on, title('Velocity vs Time'),legend,xlabel('Time [s]'), ylabel('Velocity [m/s]')
