function n = dens_profile(R,par)
    n = par(1)*exp(-((R - par(2))/par(3)).^2)+par(4)*exp(-((R-par(5))/par(6)).^2);
            end
