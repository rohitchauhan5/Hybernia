%Solves the ODE given by Yablonskiy's eqn 5. Returns time in min and
%temperature in degrees Celcius. aka temp washout for total brain 

tspan = [0 60]; %time range in minutes
mass = 100; %grams
for T_arterial = 32:37 %run through different arterial blood temperatures %24-37
    [t,heat] = ode45(@(t,temp) eqn5(t,temp,T_arterial,mass),tspan, 37*.37*mass/100); %initial heat content based on mc*delT
    figure(1)
    plot(t,heat);
    hold on; 
end
%hold off
%[t,temp]
title(sprintf('Heat Washout Curves for Cold Infusion in Healthy Brain - %d g - kJ vs time', mass));
ylabel('Heat Content of Brain (kJ)');
xlabel('time (min)');
legend('Perfusate temperature = 32','Perfusate temperature = 33', 'Perfusate temperature = 34', 'Perfusate temperature = 35', 'Perfusate temperature = 36', 'Perfusate temperature = 37');
