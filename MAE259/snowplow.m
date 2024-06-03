clear all
close all

p.Lg=35e-9;     %H
p.Rg=0.3;       %Ohm
p.l=0.04;       %m
p.c=2.998e8;    %m/s
p.r0=0.02;      %m
p.mu0=4*pi*1.0e-7;
p.M0=3e-7;      %kg
p.config=0;


s.h=1e-9;
s.T=140e-8;
s.CR=10;

t0=0;
x0=[0,p.r0,0,3e-12];   %[I,rp,up,M]

%figure(1), plot(t0,x0(2),'DisplayName','Radius'), ylim([0,0.025]), hold on, title('Radius vs Time'), xlabel('Time [s]'), ylabel('Radius [m]')
%figure(2), plot(t0,0), ylim([0,5e6],'DisplayName','Radius'), hold on, title('Voltage/Current vs Time'),xlabel('Time [s]'), ylabel('Radius [m]')


[x,t]=feval(@snowplow_RK4,@snowplow_RHS,x0,0,s,p,@snowplow_plotter);

figure(1), plot(t(:),x(2,:),'DisplayName','Radius'), ylim([0,0.025]), hold on, title('Radius vs Time'),legend,xlabel('Time [s]'), ylabel('Radius [m]')
figure(2), plot(t(:),x(1,:),'DisplayName','Current'), ylim([0,4e6]), hold on, title('Voltage/Current vs Time'),legend,xlabel('Time [s]'), ylabel('')
figure(2), plot(t(:),snowplow_Vg(t(:)),'DisplayName','Voltage'),legend,
figure(3), plot(t(:),x(3,:),'DisplayName','Velocity'), ylim([-8e5,10]), hold on, title('Velocity vs Time'),legend,xlabel('Time [s]'), ylabel('Velocity [m/s]')

%t_stop=t(end);
%x(2,end);
%[t,x] = ode45(@(t,x) snowplow_RHS(x,t,p).',[t0 t_stop],x0);


%figure(4), plot(t,x(:,2),'DisplayName','Radius'), ylim([0,0.025]), hold on, title('ode45 Radius vs Time'),legend,xlabel('Time [s]'), ylabel('Radius [m]')
%figure(5), plot(t,x(:,1),'DisplayName','Current'), ylim([0,4e6]), hold on, title('ode45 Voltage/Current vs Time'),legend,xlabel('Time [s]'), ylabel('')
%figure(5), plot(t(:),snowplow_Vg(t(:)),'DisplayName','Voltage'),legend,
%figure(6), plot(t,x(:,3),'DisplayName','Velocity'), ylim([-8e5,10]), hold on, title('ode45 Velocity vs Time'),legend,xlabel('Time [s]'), ylabel('Velocity [m/s]')


density_grid=0:1e-5:p.r0;
length(density_grid)
for i=1:length(density_grid), density(i)=snowplow_density(density_grid(i),p); end
figure(7), plot(density_grid,density), ylim([0 0.15])
