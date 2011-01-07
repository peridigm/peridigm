function [t,u] = TwoPointImplicitAnalyticalSolution(numPoints,omega)

period = 2*pi/omega;
t=0:period/numPoints:period;
u=-sin(omega*t)/omega;

end

