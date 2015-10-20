function [yout] = ForwEuler(neq, t0, tend, y0, h, func)
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
        
        ylast = ylast + h * func(tstep, ylast);
        
        ind= ind + 1;
        tstep = tstep + h;
        
        yout((ind-1) * neq + 1:ind *neq) = ylast;
    end
end