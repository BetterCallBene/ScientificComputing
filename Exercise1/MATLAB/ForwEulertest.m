neq = 3;
t0 = 0;
tend = 12;
h = 0.2;
y0 = [0, 1, 1]';
y = RungeKutta(neq, t0, tend, y0,h, @rigid);