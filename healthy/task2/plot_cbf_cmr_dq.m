tspan = [0 60]; %time range in minutes
mass = input('What is the mass in grams of healthy tissue? ')
% auc_cmr=zeros(6,2);
% auc_cbf=zeros(6,2);
% auc_dq=zeros(6,2);
% ttp_cmr=zeros(6,2); %time to plateau
% ttp_cbf=zeros(6,2);
% ttp_dq=zeros(6,2);
% j=1;

H0= 470; %  kJ/mol O2; from Yablonskiy paper,2000
Hb= 28; % kJ/mol O2; from Yablonskiy paper
p_blood = 1; % g/ml; assumed to be same as for water, Yablonskiy
c_blood = 4.178*10^-3; % specific heat in kJ/g/(degree celcius change); assumed to be same as for water, Yablonskiy

for T_arterial = 32.3:38.3 %run through different arterial blood temperatures
    [t,temp] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,38.3);
    cbf_w_time = zeros(length(t), 1);
    cmr_w_time = zeros(length(t), 1);
    dq_w_time = zeros(length(t), 1);

    for i = 1:length(temp) %through the temps as they change with time
        [cbf,cmr] = CBF_CMR_calculator(temp(i,1));
        cbf_w_time(i,1) = cbf*(mass/100); %convert cbf from ml/100g/min to ml/min
        cmr_w_time(i,1) = cmr*10^6*(mass/100); %convert cmr from mol02/100g/min to umol/min
        dq_w_time(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temp(i,1)-T_arterial))*(mass/100); %kJ/min
    end
    
    figure(2)
    plot(t,cbf_w_time);
    hold on; 
    
    figure(3)
    plot(t,cmr_w_time);
    hold on; 
    
    figure(4)
    plot(t,dq_w_time);
    hold on; 
    
%     % find integrals up to plateau
%     auc_cmr(j,1)= T_arterial;
%     value1=find_plateau(cbf_w_time, t);
%     value2=find_plateau(dq_w_time, t);
%     value=find_plateau(cmr_w_time, t);
%     auc_cmr(j,2)=trapz(t(1:value, 1),cmr_w_time(1:value,1));
%     ttp_cmr(j,1) = T_arterial;
%     ttp_cmr(j,2) = t(value);
%     auc_cbf(j,1)=T_arterial;
%     auc_cbf(j,2)=trapz(t(1:value1,1),cbf_w_time(1:value1,1));
%     ttp_cbf(j,1) = T_arterial;
%     ttp_cbf(j,2) = t(value1);
%     auc_dq(j,1)=T_arterial;
%     auc_dq(j,2)=trapz(t(1:value2,1),dq_w_time(1:value2,1));
%     ttp_dq(j,1) = T_arterial;
%     ttp_dq(j,2) = t(value2);
    j=j+1;
end
figure(2)
title('CBF vs Time for 4g Healthy Tissue with Perfusate Temps 32.3-38.3 deg C');
ylabel('CBF (ml/min)');
xlabel('time (min)');
legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3', 'Perfusate temperature = 38.3');

figure(3)
title('CMR vs Time for 4g Healthy Tissue with Perfusate Temps 32.3-38.3 deg C');
ylabel('CMR (umol02/min)');
xlabel('min');
legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3','Perfusate temperature = 38.3');

figure(4)
title('Heat Washout (dQ) vs Time for 4g Healthy Tissue with Perfusate Temps 32.3-38.3 deg C');
ylabel('Rate of Change of Heat (kJ/min)');
xlabel('min');
legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3','Perfusate temperature = 38.3');

% disp('auc cmr');
% disp(auc_cmr);
% disp('ttp cmr');
% disp(ttp_cmr);
% disp('auc cbf');
% disp(auc_cbf);
% disp('ttp cbf');
% disp(ttp_cbf);
% disp('auc dq');
% disp(auc_dq);
% disp('ttp dq');
% disp(ttp_dq);