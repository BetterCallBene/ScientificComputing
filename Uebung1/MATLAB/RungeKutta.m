function [yout] = RungeKutta(neq, t0, tend, y0, h, func)

    tstep = t0;
    nelm = ceil((tend - t0)/h) + 1;
    
    yout = zeros(nelm * neq, 1);
    done = false;
    ind = 1;
    ylast = y0;
    yout((ind-1) * neq + 1:ind *neq) = ylast;
    eps = 1e-10;
    
    while(~done)
        
        if (abs(tstep + h - tend) < eps )
            done = true;
            h = tend - tstep;
        end
        k1 = h * func(tstep, ylast);
        k2 = h * func(tstep + 1/2 * h, ylast + 1/2 * k1);
        k3 = h * func(tstep + 1/2 * h, ylast + 1/2 * k2);
        k4 = h * func(tstep + h, ylast + k3);
        
        ylast = ylast + 1/6*(k1 + 2 * k2 + 2*k3 + k4);
%       
        ind= ind + 1;
        tstep = tstep + h;
        
        yout((ind-1) * neq + 1:ind *neq) = ylast;
    end

end