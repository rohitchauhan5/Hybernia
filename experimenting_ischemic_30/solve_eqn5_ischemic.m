%Solves the ODE given by Yablonskiy's eqn 5. Returns time in min and
%temperature in degrees Celcius.

tspan = [0 60]; %time range in minutes
for T_arterial = 32.3:37.3 %run through different arterial blood temperatures
    [t,temp] = ode45(@(t,temp) eqn5_ischemic(t,temp,T_arterial),tspan, 37);
    figure(1)
    plot(t,temp);
    hold on;
end
%hold off
%[t,temp]
title('Temperature Washout Curves for Cold Infusion in Ischemic Brain');
ylabel('temperature (degrees Celcius)');
xlabel('time (min)');
legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3');