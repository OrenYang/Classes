function Burger_KS_IMEXRKCB3C_PS_FS
rng("default")
%%%%%%%%%%%%%%%%%%%% Initialize the simulation paramters (user input) %%%%%%%%%%%%%%%%%%%%
L=400; Tmax=100; N=2048; dt=0.05; PlotInt=10; alpha=1; % alpha=0 for Burgers, alpha=1 for KS
dx=L/N; x=(0:N-1)'*dx; u=0.15*randn(N,1); uhat=RC_RFFT(u,N);

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
aex43 = 1660544566939/2334033219546;
aim31 = b1; aim41 = b1; aim42=b2;  aim43=b3; aim44=b4;
aex31=b1; aex41=b1; aex42=b2;

b = [b1 b2 b3 b4]*dt;
aim = [0 0 0 0; aim21 aim22 0 0; aim31 aim32 aim33 0; aim41 aim42 aim43 aim44]*dt;
aex = [0 0 0 0; aex21 0 0 0; aex31 aex32 0 0; aex41 aex42 aex43 0]*dt;

tic
for k=1:Tmax/dt
  for rk=1:3 %%%%%%%%%%%%%%%%%%%%% ALL RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rhat = uhat;
    if (rk~=1)
      for j=1:(rk-1)
        rhat = rhat+aim(rk,j).*f(:,j)+aex(rk,j).*g(:,j);
      end
    end
    f(:,rk)=Aop.*rhat./(1-aim(rk,rk)*Aop);
    rhat = rhat+aim(rk,rk).*f(:,rk);
    rhat(fix(N/3)+1:end) = 0;   %Dealias
    r = RC_RFFTinv(rhat,N); r = -1/2*r.*r; rhat = i*kx.*RC_RFFT(r,N);
    g(:,rk)=rhat;
  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:3
    uhat = uhat+b(j)*f(:,j)+b(j)*g(:,j); % UPDATE UHAT
  end
  rs(k,:)=RC_RFFTinv(uhat,N)'; ts(k)=k*dt; % These variables are just used for plotting...
  if (mod(k,PlotInt)==0)
    pause(0.001); RC_PlotXY(x,rs(k,:),k*dt,0,L,-3,3);
    % Uncomment the lines below to make some additional interesting plots.
    % figure(2); semilogy(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([0 3 1e-8 1e-1])
    % figure(3); loglog(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([3e-2 4 1e-8 1e-1])
  end
end
timeVal = toc
timeVal/(Tmax/dt)

figure(4); rs(:,N+1)=rs(:,1); xs=[0:N]*L/N;
contour(xs,ts,rs,[.25 .75 1.25],'r-'); hold on; contour(xs,ts,rs,[-.25 -.75 -1.25],'b-.')

end % function Burger_KS_IMEXRKCB3C_PS_FS
