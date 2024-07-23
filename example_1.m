clear ; clc
% simulate the cell at a fixed voltage for a specified duration and plot
% the decay in current output caused by the degradation mechanism.

mu = 1e-2;
degradation_function = @(t) 1 + mu*t;
duration = 10*60; 

params = paper_params(degradation_function,duration,'fixed_voltage','A');
params.Verbose = true;
sol = numericalsolver(params);

%%

figure(1)
plot(sol.time/60,sol.J)

xlabel('time [minutes]')
ylabel('current density [mAcm^-^2]')
title(['current decay, \mu = ' num2str(mu,5) ' s^-^1'])


