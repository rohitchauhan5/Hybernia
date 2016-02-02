function eq_CMRO2 = CMRO2()
%----------------White matter---------------------
q0_White=4175 % q0=4175 ml/min/100g
a2=3
beta=0.1
T2=linspace(30,37) 
eq_CMRO2_White=q0_White*power(a2,(beta*(T2-37)));
figure(1);
subplot(1,2,1)
plot(T2,eq_CMRO2_White);
title('CMRO2 vs Temperature, White Matter');
xlabel('Temperature');
ylabel('CMRO2');
%----------------Grey matter----------------------
q0_Grey=16700 %q0=16700 ml/min/100g

eq_CMRO2_Grey=q0_Grey*power(a2,(beta*(T2-37)));
figure(1);
subplot(1,2,2)
plot(T2,eq_CMRO2_Grey);
title('CMRO2 vs Temperature, Grey Matter');
xlabel('Temperature');
ylabel('CMRO2');

end