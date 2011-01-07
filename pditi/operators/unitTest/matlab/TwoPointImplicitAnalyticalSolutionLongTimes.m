function [t,u] = TwoPointImplicitAnalyticalSolutionLongTimes(n,N,omega)
%%
% n -- number of periods to plot
% N -- number of points per period to plot
% omega -- natural frequency
period = 2*pi/omega;
t=0:period/N:n*period;
u=-sin(omega*t)/omega;

end