Np=201; V=[1 1.0000000001]; close all;
for k = 1:1
  switch k
    case {1}, B=[-2.55; 2.05; -2.25; 2.25];
  end
  LR(:,1)=[B(1):(B(2)-B(1))/(Np-1):B(2)]';LI(:,1)=[B(3):(B(4)-B(3))/(Np-1):B(4)]';
  for j=1:Np, for i=1:Np, L=LR(j)+sqrt(-1)*LI(i);
    switch k
      case 1, sig(i,j)=abs(1/2*(L^2-L*sqrt(-L^2-4)+2)); %SI1
    end
  end, end
  figure(k), contourf(LR,LI,1./sig,V,'k-'), colormap autumn, axis('square'), hold on
  plot([B(1) B(2)], [0,0], 'k-'), plot([0,0], [B(3) B(4)], 'k-'),
  set(gca,'FontSize',15),
  hold off;
end
