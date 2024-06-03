%% clear old data and set work directory
clearvars;
close all;
script_mpath = mfilename('fullpath');
[script_path, scipt_name] = fileparts(script_mpath);
path = script_path;
addpath(genpath(path));
cd (path);

%% ALL UNITS IN SI (MKS) except where noted

%% physical constants
c = 299792458;
mp = 1.0e-3/6.02214086e23;
mu0 = 4*pi*1.0e-7;
Ev_J = 1.602e-19;

%% circuit parameters
Cbr = 20*(58e-9/2.0);
Rbr = 15e-3;   % bricks resistance (parallel)
Lbr = 9.0e-9;  % bricks inductance (parallel)
Lpf = 10.0e-9; % feed inductance
Rpf = 6e-3;    % feed resistance
Rsh = 1.9;     % shunt (core) resistance
Vch = 110e3;   % charge voltage
Estore = 0.5*Cbr*Vch^2;

Rret = 6.0e-2; % return current radius
Hgap = 1.4e-2; % AK gap

%% load parameters
mu = 40.0; % atomic mass
n0 = 0.0e23*(4.0/mu); % peak density - CENTER
R0 = 0.0e-2; % initial radius - CENTER
DR0 = 4.0e-3; % initial width - CENTER
n1 = 3.6e22; % peak density - SHELL
R1 = 1.27e-2; % initial radius - SHELL
DR1 = 3.5e-3; % initial width - SHELL
gamma = 5/3; % polytropic index
P_coeff = 1.0e4; % pressure multiplier
Rmin = 1.0e-3; % approx. min radius
Rf = R1 + 4*DR1; %initial radius
%Rf = 15.0e-3; %initial radius
Lstart = (mu0/(2*pi))*Hgap*log(Rret/Rf); %initial load inductance

%% initial density profile
dens_params = [n0, R0, DR0, n1, R1, DR1]; % parameters to pass to the density profile function
Rgrid = 0:1e-8:Rf; % grid to calculate mass
n = dens_profile(Rgrid,dens_params); % initial denstity profile
rho = n*mu*mp; % initial mass profile
M_L = cumtrapz(Rgrid,2*pi.*Rgrid.*rho); % mass/length profile
ml = trapz(Rgrid,2*pi.*Rgrid.*rho); % total mass/length
w_bd = 0.5e-3; % initial sheath width
R_bd = Rf:1.0e-6:Rf+w_bd; % grid to calculate initial mass
N_bd = dens_profile(R_bd,dens_params); % initial sheath density
M0 = trapz(R_bd,2*pi*mu*mp.*R_bd.*N_bd); % initial sheath mass

%% simulation domain, parameters and initial conditions
t0 = 0.0; % initial time
tf = 0.8*pi*sqrt((Lbr+Lpf)*Cbr); %end time
t_span = [t0 tf]; % time limits
dt = 1.0e-11; % time step
time = t0:dt:tf; % time grid
circuit_param = [Lbr, Rbr, Cbr, Lpf, Rpf, Rsh]; % circuit parameters to pass to solver
gas_param = [mu, M0+ml, gamma]; % gas parameters to pass to solver
load_param = [Rret Hgap]; % load parameters to pass to solver
sim_param = [circuit_param gas_param load_param dt P_coeff Rf Rmin]; % simulation parameters
x0 = [Lstart -Vch*Cbr 0.0 M0 Rf 0.0 0.0]; % initial conditions
% x = [L Q I M R V Is] are the variables

x0

%% solve ODE and export variables
[t,x] = ode45(@(t,x) eqn_RLCcore_Zpinch_snowplow(t,x,sim_param,@dens_profile,dens_params),time,x0);

time_out = t0:1.0e-10:tf; % output time undersampling
L = interp1(t,x(:,1),time_out); % load inductance
Q = interp1(t,x(:,2),time_out); % charge on capacitor
I = interp1(t,x(:,3),time_out); % load current
M = interp1(t,x(:,4),time_out); % piston mass
R = interp1(t,x(:,5),time_out); % piston radius
V = interp1(t,x(:,6),time_out); % piston velocity
Is = interp1(t,x(:,7),time_out); % shunt current

%% figure 0 with initial density profile
fig0 = figure('Name','Initial density profile','units','normalized','Position',[0.02 0.05 0.95 0.85]);
set(groot,'CurrentFigure',fig0);
set(gca,'FontSize',28);
hold on;
%box on;
yyaxis left;
plot(Rgrid*1e2,rho*1e-3,'-b','LineWidth',2,'DisplayName','Initial Density');
xlim(1e2*[0 Rf]);
ylim(1e-3*[0 1.1*max(rho)]);
xlabel('Radius (cm)');
ylabel('Mass Density (g/cm^{3})');
yyaxis right;
plot(Rgrid*1e2,M_L*10,'--r','LineWidth',2,'DisplayName','Mass/Length');
xlim(1e2*[0 Rf]);
ylim(10*[0 1.1*ml]);
xlabel('Radius (cm)');
ylabel('Mass per length (g/cm)');
legend('Location','east');
hold off;

%% figure 1 with current and radius
fig1 = figure('Name','Current and Radius','units','normalized','Position',[0.02 0.05 0.95 0.85]);
set(groot,'CurrentFigure',fig1);
set(gca,'FontSize',28);
hold on;
yyaxis left;
plot(time_out*1e9,I*1e-3,'-k','LineWidth',2,'DisplayName','Current');
plot(time_out*1e9,Is*1e-3,'--k','LineWidth',0.5,'DisplayName','Shunt Current');
xlim(1e9*[time_out(1) time_out(end)]);
ylim(1.1e-3*[min(I) max(I)]);
xlabel('Time (ns)');
ylabel('Current (kA)');
yyaxis right;
plot(time_out*1e9,1e2*R,'-r','LineWidth',2,'DisplayName','Plasma radius');
ylim(1e2*[0.0 Rf]);
ylabel('Radius (cm)');
legend('Location','northwest');
hold off;
savefig(fig1,'Zpinch_I-R.fig');

%% figure 2 with radius and velocity
fig2 = figure('Name','Radius and Velocity','units','normalized','Position',[0.02 0.05 0.95 0.85]);
set(groot,'CurrentFigure',fig2);
set(gca,'FontSize',28);
hold on;
yyaxis left;
plot(time_out,R,'-k','LineWidth',2,'DisplayName','Piston Radius');
xlim([time_out(1) time_out(end)]);
ylim([0 Rf]);
xlabel('Time (s)');
ylabel('Radius (m)');
yyaxis right;
plot(time_out,V,'-r','LineWidth',2,'DisplayName','Piston Velocity');
ylim([-4 2]*1e5);
ylabel('Velocity (m/s)');
legend('Location','northwest');
hold off;
savefig(fig2,'R-V.fig');
