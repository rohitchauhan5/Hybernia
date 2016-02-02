function eq_CBF = CBF1()
a=3
beta1=0.1
T1=linspace(30,37) 
%----------------White matter---------------------
w0_White=80; %w0=80 ml/min/100gfigure(1);
eq_CBF_White = w0_White*power(a,(beta1*(T1-37)));
figure(1);
subplot(1,2,1)
plot(T1,eq_CBF_White);
title('CBF vs Temperature, White Matter');
xlabel('Temperature');
ylabel('CMRO2');
%----------------Grey matter----------------------
w0_Grey=20; %w0=20 ml/min/100g
eq_CBF_Grey = w0_Grey*power(a,(beta1*(T1-37)));

figure(2);
subplot(1,2,2)
plot(T1,eq_CBF_Grey);
title('CBF vs Temperature, Grey Matter');
xlabel('Temperature');
ylabel('CBF');
end
