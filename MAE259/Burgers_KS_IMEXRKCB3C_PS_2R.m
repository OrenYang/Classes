function Burger_KS_IMEXRKCB3C_PS_2R

%%%%%%%%%%%%%%%%%%%% Initialize the simulation paramters (user input) %%%%%%%%%%%%%%%%%%%%
L=200; Tmax=100; N=1024; dt=0.05; PlotInt=10; alpha=1; % alpha=0 for Burgers, alpha=1 for KS
dx=L/N; x=(0:N-1)'*dx; u=0.15*sin(x); uhat=RC_RFFT(u,N);
%u=0.15*randn(N,1);
%%%%%%%%%%%% Precalculate the time-stepping coefficients used in the simulation %%%%%%%%%%
kx=(2*pi/L)*[0:N/2-1]';
if alpha==0; Aop=-kx.^2; else Aop=kx.^2-kx.^4; end;
aim21 = 0;
aim22 = 3375509829940/4525919076317;
aim32 = -11712383888607531889907/32694570495602105556248;
aim33 = 566138307881/912153721139;
b1 = 0;
b2 = 673488652607/2334033219546;
b3 = 493801219040/853653026979;
b4 = 184814777513/1389668723319;
aex21 = 3375509829940/4525919076317;
aex32 = 272778623835/1039454778728;
aim43 = 1660544566939/2334033219546;

abim = [0 0 aim32-b2];
abex = [0 aex21-b1 aex32-b2];
b = [b1 b2 b3 b4];
aim = [0 aim22 aim33];

tic
for k=1:Tmax/dt
  for rk=1:3 %%%%%%%%%%%%%%%%%%%%% ALL RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (rk == 1)
      rhat = uhat;
    else
      rhat = uhat+abim(rk)*dt*Aop.*rhat+abex(rk)*dt*g(rhat,kx,N);
    end
    rhat = rhat./(1-aim(rk)*dt*Aop);
    uhat = uhat+b(rk)*dt*Aop.*rhat+b(rk)*dt.*g(rhat,kx,N);
  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rs(k,:)=RC_RFFTinv(uhat,N)'; ts(k)=k*dt; % These variables are just used for plotting...
  if (mod(k,PlotInt)==0)
    pause(0.001); RC_PlotXY(x,rs(k,:),k*dt,0,L,-1.5,1.5);
    % Uncomment the lines below to make some additional interesting plots.
    % figure(2); semilogy(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([0 3 1e-8 1e-1])
    % figure(3); loglog(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([3e-2 4 1e-8 1e-1])
  end
end
timeVal = toc
timeVal/(Tmax/dt)

figure(4); rs(:,N+1)=rs(:,1); xs=[0:N]*L/N;
contour(xs,ts,rs,[.25 .75 1.25],'r-'); hold on; contour(xs,ts,rs,[-.25 -.75 -1.25],'b-.')

end % function Burger_KS_IMEXRKCB3C_PS_2R

function [rhat] = g(rhat,kx,N)
  rhat(fix(N/3)+1:end) = 0;   %Dealias
  r = RC_RFFTinv(rhat,N); r = -1/2*r.*r; rhat = i*kx.*RC_RFFT(r,N);
end
