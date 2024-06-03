function snowplowPlotter(xo,xn,to,tn, ooo, xxx)
figure(1), plot([to tn],[xo(2) xn(2)]), hold on, pause(0.00000001)
figure(2), plot([to tn],[ooo xxx]), hold on, pause(0.00000001)
figure(2), plot([to tn],[xo(1) xn(1)]), hold on, pause(0.00000001)
end
