% T_prompt= 'Enter the arterial temperature'
% T_Arterial=input(T_prompt);
% C_prompt= ' Enter the heat capacity of the tissue'
% C_Tissue=(C_prompt);
% % C_Tissue is heat capacity of the tissue, which is generally equal to 3630
% % J/kg/deg C
% rho_prompt= ' Enter the density of the tissue'
% rho=input(rho_prompt);
% Cb_prompt= ' Enter the Cb of the tissue'
% Cb=input(Cb_prompt);
% w0_prompt= ' Enter the w0 of the tissue';
% w0=input(w0_prompt);
% q0_prompt= ' Enter the q0 of the tissue';
% q0=input(q0_prompt);
% beta1_prompt= ' Enter the beta constant of CBF';
% beta1=input(beta1_prompt);
% beta_prompt= ' Enter the beta constant of CMRO2';
% beta=input(beta_prompt);
% a_prompt= ' Enter the alpha constant of CBF';
% a=input(a_prompt);
% a2_prompt= ' Enter the alpha constant of CMRO2';
% a2=input(a2_prompt);

delta_H0=470 % delta_H0 is 470kJ per mol O2 from Yablonskiy paper,2000
delta_Hb=28 % delta_Hb is 28kJ per mol O2 from Yablonskiy paper,2000
a2=3;%unitless
beta=0.1;%unitless
rho=1030; % kg/m3
Cb=420; %J
a=3;%unitless
beta1=0.1;%unitless
T_Arterial=37;%Unit-Degree Celcius
C_Tissue=3700;%J/kg/degree Celcius
T=linspace(30,37);% Unit-Degree Celcius
% --------------------------- White Matter------------------------------
q0_White=4175;% Unit- W/m3
w0_White=20;%ml/min/100g
for i = 1:100
T(i)
T_dot_White(i)=((delta_H0-delta_Hb)*q0_White*power(a2,(beta*(T(i)-37)))-rho*Cb*w0_White*power(a,(beta1*(T(i)-37)))*(T(i)-T_Arterial))/C_Tissue;

end;
figure(1)
subplot(1,2,1)
%make here your first plot
plot(T,T_dot_White);
title('Plot of Rate of Temperature change with Temperature for white matter');
xlabel('Temperature');
ylabel('Rate of Temperature change');

% --------------------------- Grey Matter------------------------------
q0_Grey=4175;% Unit- W/m3
a2=3;%unitless
w0_Grey=20;%ml/min/100g

for i = 1:100
T(i)
T_dot_Grey(i)=((delta_H0-delta_Hb)*q0_Grey*power(a2,(beta*(T(i)-37)))-rho*Cb*w0_Grey*power(a,(beta1*(T(i)-37)))*(T(i)-T_Arterial))/C_Tissue;

end;
figure(2)
subplot(1,2,2)
%make here your second plot
plot(T,T_dot_Grey);
title('Plot of Rate of Temperature change with Temperature for grey matter');
xlabel('Temperature');
ylabel('Rate of Temperature change');

