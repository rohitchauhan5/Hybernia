tspan = [0 60]; %time range in minutes
mass = input('What is the mass in grams of healthy tissue?')
auc_cmr=zeros(5,2);
auc_cbf=zeros(5,2);
j=1;
for T_arterial = 32.3:36.3 %run through different arterial blood temperatures
    [t,temp] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,37);
    cbf_w_time = zeros(length(t), 1);
    cmr_w_time = zeros(length(t), 1);
    %total_brain_bf_w_time = zeros(length(t),1);
    for i = 1:length(temp) %through the temps as they change with time
        [cbf,cmr] = CBF_CMR_calculator(temp(i,1));
        cbf_w_time(i,1) = cbf*(mass/100); %convert cbf from ml/100g/min to ml/min
        cmr_w_time(i,1) = cmr*10^6*(mass/100); %convert cmr from mol02/100g/min to umol/min
        %total_brain_bf_w_time(i,1) = cbf*5; %assuming brain of 500g, simply cbf*5
    end
    
    figure(2)
    plot(t,cbf_w_time);
    hold on; 
    
    figure(3)
    plot(t,cmr_w_time);
    hold on; 
    auc_cmr(j,1)=T_arterial;
    value=find_plateau(cmr_w_time);
    value1=find_plateau(cbf_w_time);
    auc_cmr(j,2)=trapz(t(1:value, 1),cmr_w_time(1:value,1));
    auc_cbf(j,1)=T_arterial;
    auc_cbf(j,2)=trapz(t(1:value1,1),cbf_w_time(1:value1,1));
    j=j+1;
%     figure(4)
%     plot(t, total_brain_bf_w_time);
%     hold on;
end
figure(2)
title('CBF vs Time for 500g Healthy Tissue with Perfusate Temps 32.3-36.3 deg C');
ylabel('CBF (ml/min)');
xlabel('time (min)');
legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3');

figure(3)
title('CMR vs Time for 500g Healthy Tissue with Perfusate Temps 32.3-36.3 deg C');
ylabel('CMR (umol02/min)');
xlabel('min');
legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3');

% figure(4)
% title('Total Brain Blood Flow vs Time for Perfusate Temps 24-37 deg C');
% ylabel('Blood Flow (ml/min)');
% xlabel('min');
% legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3');
disp(auc_cmr);
disp(auc_cbf);

