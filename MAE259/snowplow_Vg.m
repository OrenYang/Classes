function f=snowplow_Vg(t)
Vg0=2.5e6;    %V
tg=95e-9;     %s
f=Vg0.*(t/tg).^2.*exp(1-(t/tg).^2);

end
