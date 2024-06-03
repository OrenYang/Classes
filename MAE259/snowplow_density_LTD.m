function n = snowplow_density_LTD(R,p)
    n = p.n0*exp(-((R - p.R0)/p.DR0).^2)+p.n1*exp(-((R-p.R1)/p.DR1).^2);
end
