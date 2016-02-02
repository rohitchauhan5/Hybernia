%   displays auc and ttp for cbf, cmr, and dq
%   do not call in a loop, will overwrite itself
%   cbf_w_time, cmr_w_time, dq overwrite with each iteration of T_arterial

base_cbf = input('What is baseline cbf? ')
base_cmr = input('What is baseline cmr? ')

tspan = [0 60]; %time range in minutes
mass = 500; %input('What is the mass in grams of healthy tissue? ')
auc_cmr=zeros(6,2);
auc_cbf=zeros(6,2);
auc_dq=zeros(6,2);
ttp_cmr=zeros(6,2); %time to plateau
ttp_cbf=zeros(6,2);
ttp_dq=zeros(6,2);
j=1;

H0= 470; %  kJ/mol O2; from Yablonskiy paper,2000
Hb= 28; % kJ/mol O2; from Yablonskiy paper
p_blood = 1; % g/ml; assumed to be same as for water, Yablonskiy
c_blood = 4.178*10^-3; % specific heat in kJ/g/(degree celcius change); assumed to be same as for water, Yablonskiy

% *******EDIT LINE BELOW TO HAVE ONLY 1 TEMP FOR T_ARTERIAL IF WANT TO VIEW
% cbf, cmr, and dq FOR THAT T_ARTERIAL********
for T_arterial = 32.3:37.3 %run through different arterial blood temperatures
    [t,temp] = ode45(@(t,temp) eqn5(t,temp,T_arterial, base_cbf, base_cmr*10^-6),tspan,37.3);
    cbf_w_time = zeros(length(t), 1);
    cmr_w_time = zeros(length(t), 1); 
    dq_w_time = zeros(length(t), 1);
    
    for i = 1:length(temp) %through the temps as they change with time
        [cbf,cmr] = CBF_CMR_calculator_give_pars(temp(i,1), base_cbf, base_cmr*10^(-6));
        cbf_w_time(i,1) = cbf*(mass/100); %convert cbf from ml/100g/min to ml/min
        cmr_w_time(i,1) = cmr*10^6*(mass/100); %convert cmr from mol02/100g/min to umol/min
        dq_w_time(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temp(i,1)-T_arterial))*(mass/100); %kJ/min
    end

% find integrals up to plateau
auc_cmr(j,1)= T_arterial;
value1=find_plateau(cbf_w_time, t);
value2=find_plateau(dq_w_time, t);
value=find_plateau(cmr_w_time, t);
auc_cmr(j,2)=trapz(t(1:value, 1),cmr_w_time(1:value,1));
ttp_cmr(j,1) = T_arterial;
ttp_cmr(j,2) = t(value);
auc_cbf(j,1)=T_arterial;
auc_cbf(j,2)=trapz(t(1:value1,1),cbf_w_time(1:value1,1));
ttp_cbf(j,1) = T_arterial;
ttp_cbf(j,2) = t(value1);
auc_dq(j,1)=T_arterial;
auc_dq(j,2)=trapz(t(1:value2,1),dq_w_time(1:value2,1));
ttp_dq(j,1) = T_arterial;
ttp_dq(j,2) = t(value2);

j=j+1;
end

disp('auc cmr');
disp(auc_cmr);
disp('ttp cmr');
disp(ttp_cmr);
disp('auc cbf');
disp(auc_cbf);
disp('ttp cbf');
disp(ttp_cbf);
disp('auc dq');
disp(auc_dq);
disp('ttp dq');
disp(ttp_dq);

disp('t')
disp(t)
disp('temp')
disp(temp)
disp('cbf')
disp(cbf_w_time)
disp('cmr')
disp(cmr_w_time)
disp('dq')
disp(dq_w_time)
end

