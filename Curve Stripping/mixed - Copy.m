tspan = [0:1:60]; %time range in minutes, specify intervals so will be same soln for healthy and ischemic eqn5
% massh = input('What is the mass in grams of healthy tissue? ')
% massi = input('What is the mass in grams of ischemic tissue? ')
% auc_cmr36=zeros(22,2);
% auc_cbf33=zeros(22,2);
% auc_dq32=zeros(22,2);
% auc_q32=zeros(22,2);
ttp_cmr=zeros(22,2); %time to plateau
ttp_cbf=zeros(22,2);
ttp_dq=zeros(22,2);
%Rohit curve stripping
value1=zeros(6,22);
value2=zeros(6,22);
value3=zeros(6,22);
value4=zeros(6,22);
% for 37.3 degree C
b1lm_cbf_w_t37=zeros(22,2);
b2lm_cbf_w_t37=zeros(22,2);
b1lm_cmr_w_t37=zeros(22,2);
b2lm_cmr_w_t37=zeros(22,2);
b1lm_dq_w_t37=zeros(22,2);
b2lm_dq_w_t37=zeros(22,2);
b1lm_q_w_t37=zeros(22,2);
b2lm_q_w_t37=zeros(22,2);
auc_cmr37=zeros(22,2);
auc_cbf37=zeros(22,2);
auc_dq37=zeros(22,2);
auc_q37=zeros(22,2);

% for 36.3 degree C
b1lm_cbf_w_t36=zeros(22,2);
b2lm_cbf_w_t36=zeros(22,2);
b1lm_cmr_w_t36=zeros(22,2);
b2lm_cmr_w_t36=zeros(22,2);
b1lm_dq_w_t36=zeros(22,2);
b2lm_dq_w_t36=zeros(22,2);
b1lm_q_w_t36=zeros(22,2);
b2lm_q_w_t36=zeros(22,2);
auc_cmr36=zeros(22,2);
auc_cbf36=zeros(22,2);
auc_dq36=zeros(22,2);
auc_q36=zeros(22,2);


% for 35.3 degree C
b1lm_cbf_w_t35=zeros(22,2);
b2lm_cbf_w_t35=zeros(22,2);
b1lm_cmr_w_t35=zeros(22,2);
b2lm_cmr_w_t35=zeros(22,2);
b1lm_dq_w_t35=zeros(22,2);
b2lm_dq_w_t35=zeros(22,2);
b1lm_q_w_t35=zeros(22,2);
b2lm_q_w_t35=zeros(22,2);
auc_cmr35=zeros(22,2);
auc_cbf35=zeros(22,2);
auc_dq35=zeros(22,2);
auc_q35=zeros(22,2);


% for 34.3 degree C
b1lm_cbf_w_t34=zeros(22,2);
b2lm_cbf_w_t34=zeros(22,2);
b1lm_cmr_w_t34=zeros(22,2);
b2lm_cmr_w_t34=zeros(22,2);
b1lm_dq_w_t34=zeros(22,2);
b2lm_dq_w_t34=zeros(22,2);
b1lm_q_w_t34=zeros(22,2);
b2lm_q_w_t34=zeros(22,2);
auc_cmr34=zeros(22,2);
auc_cbf34=zeros(22,2);
auc_dq34=zeros(22,2);
auc_q34=zeros(22,2);

% for 33.3 degree C
b1lm_cbf_w_t33=zeros(22,2);
b2lm_cbf_w_t33=zeros(22,2);
b1lm_cmr_w_t33=zeros(22,2);
b2lm_cmr_w_t33=zeros(22,2);
b1lm_dq_w_t33=zeros(22,2);
b2lm_dq_w_t33=zeros(22,2);
b1lm_q_w_t33=zeros(22,2);
b2lm_q_w_t33=zeros(22,2);
auc_cmr33=zeros(22,2);
auc_cbf33=zeros(22,2);
auc_dq33=zeros(22,2);
auc_q33=zeros(22,2);

% % for 32.3 degree C
b1lm_cbf_w_t32=zeros(22,2);
b2lm_cbf_w_t32=zeros(22,2);
b1lm_cmr_w_t32=zeros(22,2);
b2lm_cmr_w_t32=zeros(22,2);
b1lm_dq_w_t32=zeros(22,2);
b2lm_dq_w_t32=zeros(22,2);
b1lm_q_w_t32=zeros(22,2);
b2lm_q_w_t32=zeros(22,2);
auc_cbf32=zeros(22,2);
auc_cmr32=zeros(22,2);
auc_dq32=zeros(22,2);
auc_q32=zeros(22,2);


H0= 470; %  kJ/mol O2; from Yablonskiy paper,2000
Hb= 28; % kJ/mol O2; from Yablonskiy paper
p_blood = 1; % g/ml; assumed to be same as for water, Yablonskiy
c_blood = 4.178*10^-3; % specific heat in kJ/g/(degree celcius change); assumed to be same as for water, Yablonskiy

massi=zeros(22,1);
massh=zeros(22,1);

for k=1:21 % for running a loop 20 times, for every 25 gm change in healthy and ischemic brain from 500g healthy to 0g healthy
    massh(k,1)=500-(k-1)*25;
    massi(k,1)=(k-1)*25;
    disp(k);
T_arterial = 32.3 %run through different arterial blood temperatures
    [t,temph] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,37);
    [t,tempi] = ode45(@(t,temp) eqn5_ischemic(t,temp,T_arterial),tspan,37);
    cbf_w_timeh = zeros(length(t), 1);
    cmr_w_timeh = zeros(length(t), 1);
    dq_w_timeh = zeros(length(t), 1);
    q_w_time = zeros(length(t),1); % added to calculate integral of dq/dt over t
    cbf_w_timei = zeros(length(t), 1);
    cmr_w_timei = zeros(length(t), 1);
    dq_w_timei = zeros(length(t), 1);

    for i = 1:length(temph) % healthy tissue 
        [cbf,cmr] = CBF_CMR_calculator(temph(i,1));
        cbf_w_timeh(i,1) = cbf*(massh(k,1)/100); %convert cbf from ml/100g/min to ml/min
        cmr_w_timeh(i,1) = cmr*10^6*(massh(k,1)/100); %convert cmr from mol02/100g/min to umol/min
        dq_w_timeh(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temph(i,1)-T_arterial))*(massh(k,1)/100); %kJ/min
    end
    
    
    for i = 1:length(tempi) % ischemic tissue
        [normal_cbf, normal_cmr] = CBF_CMR_calculator(tempi(i,1));
        % because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
        % per 100 g if ASSUME density of brain = p_water = 1 g/ml)
        if(normal_cbf >= 5)
            cbf = 5;
            cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
%             temp=
        else
            cbf = normal_cbf;
            %should w0 change? no
            cmr = normal_cmr;
        end
        cbf_w_timei(i,1) = cbf*(massi(k,1)/100); %convert cbf from ml/100g/min to ml/min
        cmr_w_timei(i,1) = cmr*10^6*(massi(k,1)/100); %convert cmr from mol02/100g/min to umol/min
        dq_w_timei(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(tempi(i,1)-T_arterial))*(massi(k,1)/100); %kJ/min
    end
    
    cbf_w_time = cbf_w_timei + cbf_w_timeh;
    cmr_w_time = cmr_w_timei + cmr_w_timeh;
    dq_w_time = dq_w_timei + dq_w_timeh;
    q_w_time = cumsum(dq_w_time);
    
    
    figure(1)
    plot(t,cbf_w_time);
    hold on; 
    
    
    figure(2)
    plot(t,cmr_w_time);
    hold on; 
    
    figure(3)
    plot(t,dq_w_time);
    hold on; 
    
    figure(4)
    plot(t,q_w_time);
    hold on
    
    value1(1,k)=find_plateau(cbf_w_time, t)
    value2(1,k)=find_plateau(cmr_w_time, t);
    value3(1,k)=find_plateau(dq_w_time, t);
    value4(1,k)=find_plateau(q_w_time, t);
    
    % Rohit work-curve stripping
    % CBF-Value1
    md1_cbf_w_t=fitlm(t(1:value1(1,k),1),cbf_w_time(1:value1(1,k),1));
    x_cbf_w_t=md1_cbf_w_t.Coefficients.Estimate;
    b1lm_cbf_w_t32(k,:)=x_cbf_w_t.';
    
    md2_cbf_w_t=fitlm(t(value1(1,k):60,1),cbf_w_time(value1(1,k):60,1));
    k_cbf_w_t=md2_cbf_w_t.Coefficients.Estimate;
    b2lm_cbf_w_t32(k,:)=k_cbf_w_t.';
    
    % CMR-Value2
    md1_cmr_w_t=fitlm(t(1:value2(1,k),1),cmr_w_time(1:value2(1,k),1));
    x_cmr_w_t=md1_cmr_w_t.Coefficients.Estimate;
    b1lm_cmr_w_t32(k,:)=x_cmr_w_t.';
    
    md2_cmr_w_t=fitlm(t(value2(1,k):60,1),cmr_w_time(value2(1,k):60,1));
    k_cmr_w_t=md2_cmr_w_t.Coefficients.Estimate;
    b2lm_cmr_w_t32(k,:)=k_cmr_w_t.';
    
    % dQ-Value3
    md1_dq_w_t=fitlm(t(1:value3(1,k),1),dq_w_time(1:value3(1,k),1));
    x_dq_w_t=md1_dq_w_t.Coefficients.Estimate;
    b1lm_dq_w_t32(k,:)=x_dq_w_t.';
    
    md2_dq_w_t=fitlm(t(value3(1,k):60,1),dq_w_time(value3(1,k):60,1));
    k_dq_w_t=md2_dq_w_t.Coefficients.Estimate;
    b2lm_dq_w_t32(k,:)=k_dq_w_t.';
    
    % Q-Value4
    md1_q_w_t=fitlm(t(1:value4(1,k),1),q_w_time(1:value4(1,k),1));
    x_q_w_t=md1_q_w_t.Coefficients.Estimate;
    b1lm_q_w_t32(k,:)=x_q_w_t.';
    
    md2_q_w_t=fitlm(t(value4(1,k):60,1),q_w_time(value4(1,k):60,1));
    k_q_w_t=md2_q_w_t.Coefficients.Estimate;
    b2lm_q_w_t32(k,:)=k_q_w_t.';
    
   
    %find auc integrals up to plateau for different curves, i.e.,
    %cbf, cmr, dqdt and q
    
    % Auc for cbf curve
    auc_cbf32(k,1)=massi(k,1);
    auc_cbf32(k,2)=trapz(t(1:value1(1,k),1),cbf_w_time(1:value1(1,k),1));
    
    % Auc for cmr curve
    auc_cmr32(k,1)=massi(k,1);
    auc_cmr32(k,2)=trapz(t(1:value2(1,k),1),cmr_w_time(1:value2(1,k),1));
    
    % Auc for dq curve
    auc_dq32(k,1)=massi(k,1);
    auc_dq32(k,2)=trapz(t(1:value3(1,k),1),cmr_w_time(1:value3(1,k),1));
    
    % Auc for q curve
    auc_q32(k,1)=massi(k,1);
    auc_q32(k,2)=trapz(t(1:value4(1,k),1),q_w_time(1:value4(1,k),1));
   
% 
% T_arterial = 33.3 %run through different arterial blood temperatures
%     [t,temph] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,37);
%     [t,tempi] = ode45(@(t,temp) eqn5_ischemic(t,temp,T_arterial),tspan,37);
%     cbf_w_timeh = zeros(length(t), 1);
%     cmr_w_timeh = zeros(length(t), 1);
%     dq_w_timeh = zeros(length(t), 1);
%     q_w_time = zeros(length(t),1); % added to calculate integral of dq/dt over t
%     cbf_w_timei = zeros(length(t), 1);
%     cmr_w_timei = zeros(length(t), 1);
%     dq_w_timei = zeros(length(t), 1);
% 
%     for i = 1:length(temph) % healthy tissue 
%         [cbf,cmr] = CBF_CMR_calculator(temph(i,1));
%         cbf_w_timeh(i,1) = cbf*(massh/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timeh(i,1) = cmr*10^6*(massh/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timeh(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temph(i,1)-T_arterial))*(massh/100); %kJ/min
%     end
%     
%     
%     for i = 1:length(tempi) % ischemic tissue
%         [normal_cbf, normal_cmr] = CBF_CMR_calculator(tempi(i,1));
%         % because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
%         % per 100 g if ASSUME density of brain = p_water = 1 g/ml)
%         if(normal_cbf >= 15)
%             cbf = 15;
%             cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
% %             temp=
%         else
%             cbf = normal_cbf;
%             %should w0 change? no
%             cmr = normal_cmr;
%         end
%         cbf_w_timei(i,1) = cbf*(massi/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timei(i,1) = cmr*10^6*(massi/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timei(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(tempi(i,1)-T_arterial))*(massi/100); %kJ/min
%     end
%     
%     cbf_w_time = cbf_w_timei + cbf_w_timeh;
%     cmr_w_time = cmr_w_timei + cmr_w_timeh;
%     dq_w_time = dq_w_timei + dq_w_timeh;
%     q_w_time = cumsum(dq_w_time);
%     
%     
%     figure(1)
%     plot(t,cbf_w_time);
%     hold on; 
%     
%     
%     figure(2)
%     plot(t,cmr_w_time);
%     hold on; 
%     
%     figure(3)
%     plot(t,dq_w_time);
%     hold on; 
%     
%     figure(4)
%     plot(t,q_w_time);
%     hold on
%     
%     value1(2,k)=find_plateau(cbf_w_time, t);
%     value2(2,k)=find_plateau(cmr_w_time, t);
%     value3(2,k)=find_plateau(dq_w_time, t);
%     value4(2,k)=find_plateau(q_w_time, t);
%     
%     % Rohit work-curve stripping
%     % CBF-Value1
%     md1_cbf_w_t=fitlm(t(1:value1(2,k),1),cbf_w_time(1:value1(2,k),1));
%     x_cbf_w_t=md1_cbf_w_t.Coefficients.Estimate;
%     b1lm_cbf_w_t33(k,:)=x_cbf_w_t.';
%     
%     md2_cbf_w_t=fitlm(t(value1(2,k):60,1),cbf_w_time(value1(2,k):60,1));
%     k_cbf_w_t=md2_cbf_w_t.Coefficients.Estimate;
%     b2lm_cbf_w_t33(k,:)=k_cbf_w_t.';
%     
%     % CMR-Value2
%     md1_cmr_w_t=fitlm(t(1:value2(2,k),1),cmr_w_time(1:value2(2,k),1));
%     x_cmr_w_t=md1_cmr_w_t.Coefficients.Estimate;
%     b1lm_cmr_w_t33(k,:)=x_cmr_w_t.';
%     
%     md2_cmr_w_t=fitlm(t(value2(2,k):60,1),cmr_w_time(value2(2,k):60,1));
%     k_cmr_w_t=md2_cmr_w_t.Coefficients.Estimate;
%     b2lm_cmr_w_t33(k,:)=k_cmr_w_t.';
%     
%     % dQ-Value3
%     md1_dq_w_t=fitlm(t(1:value3(2,k),1),dq_w_time(1:value3(2,k),1));
%     x_dq_w_t=md1_dq_w_t.Coefficients.Estimate;
%     b1lm_dq_w_t33(k,:)=x_dq_w_t.';
%     
%     md2_dq_w_t=fitlm(t(value3(2,k):60,1),dq_w_time(value3(2,k):60,1));
%     k_dq_w_t=md2_dq_w_t.Coefficients.Estimate;
%     b2lm_dq_w_t33(k,:)=k_dq_w_t.';
%     
%     % Q-Value4
%     md1_q_w_t=fitlm(t(1:value4(2,k),1),q_w_time(1:value4(2,k),1));
%     x_q_w_t=md1_q_w_t.Coefficients.Estimate;
%     b1lm_q_w_t33(k,:)=x_q_w_t.';
%     
%     md2_q_w_t=fitlm(t(value4(2,k):60,1),q_w_time(value4(2,k):60,1));
%     k_q_w_t=md2_q_w_t.Coefficients.Estimate;
%     b2lm_q_w_t33(k,:)=k_q_w_t.';
%     
%    
%     %find auc integrals up to plateau for different curves, i.e.,
%     %cbf, cmr, dqdt and q
%     
%     % Auc for cbf curve
%     auc_cbf33(k,1)=massi;
%     auc_cbf33(k,2)=trapz(t(1:value1(2,k),1),cbf_w_time(1:value1(2,k),1));
%     
%     % Auc for cmr curve
%     auc_cmr33(k,1)=massi;
%     auc_cmr33(k,2)=trapz(t(1:value2(2,k),1),cmr_w_time(1:value2(2,k),1));
%     
%     % Auc for dq curve
%     auc_dq33(k,1)=massi;
%     auc_dq33(k,2)=trapz(t(1:value3(2,k),1),cmr_w_time(1:value3(2,k),1));
%     
%     % Auc for q curve
%     auc_q33(k,1)=massi;
%     auc_q33(k,2)=trapz(t(1:value4(2,k),1),q_w_time(1:value4(2,k),1));
%    
% 
%     
%    
%   
% T_arterial = 34.3 %run through different arterial blood temperatures
%     [t,temph] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,37);
%     [t,tempi] = ode45(@(t,temp) eqn5_ischemic(t,temp,T_arterial),tspan,37);
%     cbf_w_timeh = zeros(length(t), 1);
%     cmr_w_timeh = zeros(length(t), 1);
%     dq_w_timeh = zeros(length(t), 1);
%     q_w_time = zeros(length(t),1); % added to calculate integral of dq/dt over t
%     cbf_w_timei = zeros(length(t), 1);
%     cmr_w_timei = zeros(length(t), 1);
%     dq_w_timei = zeros(length(t), 1);
% 
%     for i = 1:length(temph) % healthy tissue 
%         [cbf,cmr] = CBF_CMR_calculator(temph(i,1));
%         cbf_w_timeh(i,1) = cbf*(massh/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timeh(i,1) = cmr*10^6*(massh/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timeh(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temph(i,1)-T_arterial))*(massh/100); %kJ/min
%     end
%     
%     
%     for i = 1:length(tempi) % ischemic tissue
%         [normal_cbf, normal_cmr] = CBF_CMR_calculator(tempi(i,1));
%         % because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
%         % per 100 g if ASSUME density of brain = p_water = 1 g/ml)
%         if(normal_cbf >= 15)
%             cbf = 15;
%             cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
% %             temp=
%         else
%             cbf = normal_cbf;
%             %should w0 change? no
%             cmr = normal_cmr;
%         end
%         cbf_w_timei(i,1) = cbf*(massi/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timei(i,1) = cmr*10^6*(massi/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timei(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(tempi(i,1)-T_arterial))*(massi/100); %kJ/min
%     end
%     
%     cbf_w_time = cbf_w_timei + cbf_w_timeh;
%     cmr_w_time = cmr_w_timei + cmr_w_timeh;
%     dq_w_time = dq_w_timei + dq_w_timeh;
%     q_w_time = cumsum(dq_w_time);
%     
%     
%     figure(1)
%     plot(t,cbf_w_time);
%     hold on; 
%     
%     
%     figure(2)
%     plot(t,cmr_w_time);
%     hold on; 
%     
%     figure(3)
%     plot(t,dq_w_time);
%     hold on; 
%     
%     figure(4)
%     plot(t,q_w_time);
%     hold on
%     
%     value1(3,k)=find_plateau(cbf_w_time, t);
%     value2(3,k)=find_plateau(cmr_w_time, t);
%     value3(3,k)=find_plateau(dq_w_time, t);
%     value4(3,k)=find_plateau(q_w_time, t);
%     
%     % Rohit work-curve stripping
%     % CBF-Value1
%     md1_cbf_w_t=fitlm(t(1:value1(3,k),1),cbf_w_time(1:value1(3,k),1));
%     x_cbf_w_t=md1_cbf_w_t.Coefficients.Estimate;
%     b1lm_cbf_w_t34(k,:)=x_cbf_w_t.';
%     
%     md2_cbf_w_t=fitlm(t(value1(3,k):60,1),cbf_w_time(value1(3,k):60,1));
%     k_cbf_w_t=md2_cbf_w_t.Coefficients.Estimate;
%     b2lm_cbf_w_t34(k,:)=k_cbf_w_t.';
%     
%     % CMR-Value2
%     md1_cmr_w_t=fitlm(t(1:value2(3,k),1),cmr_w_time(1:value2(3,k),1));
%     x_cmr_w_t=md1_cmr_w_t.Coefficients.Estimate;
%     b1lm_cmr_w_t34(k,:)=x_cmr_w_t.';
%     
%     md2_cmr_w_t=fitlm(t(value2(3,k):60,1),cmr_w_time(value2(3,k):60,1));
%     k_cmr_w_t=md2_cmr_w_t.Coefficients.Estimate;
%     b2lm_cmr_w_t34(k,:)=k_cmr_w_t.';
%     
%     % dQ-Value3
%     md1_dq_w_t=fitlm(t(1:value3(3,k),1),dq_w_time(1:value3(3,k),1));
%     x_dq_w_t=md1_dq_w_t.Coefficients.Estimate;
%     b1lm_dq_w_t34(k,:)=x_dq_w_t.';
%     
%     md2_dq_w_t=fitlm(t(value3(3,k):60,1),dq_w_time(value3(3,k):60,1));
%     k_dq_w_t=md2_dq_w_t.Coefficients.Estimate;
%     b2lm_dq_w_t34(k,:)=k_dq_w_t.';
%     
%     % Q-Value4
%     md1_q_w_t=fitlm(t(1:value4(3,k),1),q_w_time(1:value4(3,k),1));
%     x_q_w_t=md1_q_w_t.Coefficients.Estimate;
%     b1lm_q_w_t34(k,:)=x_q_w_t.';
%     
%     md2_q_w_t=fitlm(t(value4(3,k):60,1),q_w_time(value4(3,k):60,1));
%     k_q_w_t=md2_q_w_t.Coefficients.Estimate;
%     b2lm_q_w_t34(k,:)=k_q_w_t.';
%     
%    
%     %find auc integrals up to plateau for different curves, i.e.,
%     %cbf, cmr, dqdt and q
%     
%     % Auc for cbf curve
%     auc_cbf34(k,1)=massi;
%     auc_cbf34(k,2)=trapz(t(1:value1(3,k),1),cbf_w_time(1:value1(3,k),1));
%     
%     % Auc for cmr curve
%     auc_cmr34(k,1)=massi;
%     auc_cmr34(k,2)=trapz(t(1:value2(3,k),1),cmr_w_time(1:value2(3,k),1));
%     
%     % Auc for dq curve
%     auc_dq34(k,1)=massi;
%     auc_dq34(k,2)=trapz(t(1:value3(3,k),1),cmr_w_time(1:value3(3,k),1));
%     
%     % Auc for q curve
%     auc_q34(k,1)=massi;
%     auc_q34(k,2)=trapz(t(1:value4(3,k),1),q_w_time(1:value4(3,k),1));
%  
%    
% 
% T_arterial = 35.3 %run through different arterial blood temperatures
%     [t,temph] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,37);
%     [t,tempi] = ode45(@(t,temp) eqn5_ischemic(t,temp,T_arterial),tspan,37);
%     cbf_w_timeh = zeros(length(t), 1);
%     cmr_w_timeh = zeros(length(t), 1);
%     dq_w_timeh = zeros(length(t), 1);
%     q_w_time = zeros(length(t),1); % added to calculate integral of dq/dt over t
%     cbf_w_timei = zeros(length(t), 1);
%     cmr_w_timei = zeros(length(t), 1);
%     dq_w_timei = zeros(length(t), 1);
% 
%     for i = 1:length(temph) % healthy tissue 
%         [cbf,cmr] = CBF_CMR_calculator(temph(i,1));
%         cbf_w_timeh(i,1) = cbf*(massh/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timeh(i,1) = cmr*10^6*(massh/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timeh(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temph(i,1)-T_arterial))*(massh/100); %kJ/min
%     end
%     
%     
%     for i = 1:length(tempi) % ischemic tissue
%         [normal_cbf, normal_cmr] = CBF_CMR_calculator(tempi(i,1));
%         % because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
%         % per 100 g if ASSUME density of brain = p_water = 1 g/ml)
%         if(normal_cbf >= 15)
%             cbf = 15;
%             cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
% %             temp=
%         else
%             cbf = normal_cbf;
%             %should w0 change? no
%             cmr = normal_cmr;
%         end
%         cbf_w_timei(i,1) = cbf*(massi/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timei(i,1) = cmr*10^6*(massi/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timei(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(tempi(i,1)-T_arterial))*(massi/100); %kJ/min
%     end
%     
%     cbf_w_time = log(cbf_w_timei + cbf_w_timeh);
%     cmr_w_time = cmr_w_timei + cmr_w_timeh;
%     dq_w_time = dq_w_timei + dq_w_timeh;
%     q_w_time = cumsum(dq_w_time);
%     
%     
%     figure(1)
%     plot(t,cbf_w_time);
%     hold on; 
%     
%     
%     figure(2)
%     plot(t,cmr_w_time);
%     hold on; 
%     
%     figure(3)
%     plot(t,dq_w_time);
%     hold on; 
%     
%     figure(4)
%     plot(t,q_w_time);
%     hold on
%     
%     value1(4,k)=find_plateau(cbf_w_time, t);
%     value2(4,k)=find_plateau(cmr_w_time, t);
%     value3(4,k)=find_plateau(dq_w_time, t);
%     value4(4,k)=find_plateau(q_w_time, t);
%     
%     % Rohit work-curve stripping
%     % CBF-Value1
%     md1_cbf_w_t=fitlm(t(1:value1(4,k),1),cbf_w_time(1:value1(4,k),1));
%     x_cbf_w_t=md1_cbf_w_t.Coefficients.Estimate;
%     b1lm_cbf_w_t35(k,:)=x_cbf_w_t.';
%     
%     md2_cbf_w_t=fitlm(t(value1(4,k):60,1),cbf_w_time(value1(4,k):60,1));
%     k_cbf_w_t=md2_cbf_w_t.Coefficients.Estimate;
%     b2lm_cbf_w_t35(k,:)=k_cbf_w_t.';
%     
%     % CMR-Value2
%     md1_cmr_w_t=fitlm(t(1:value2(4,k),1),cmr_w_time(1:value2(4,k),1));
%     x_cmr_w_t=md1_cmr_w_t.Coefficients.Estimate;
%     b1lm_cmr_w_t35(k,:)=x_cmr_w_t.';
%     
%     md2_cmr_w_t=fitlm(t(value2(4,k):60,1),cmr_w_time(value2(4,k):60,1));
%     k_cmr_w_t=md2_cmr_w_t.Coefficients.Estimate;
%     b2lm_cmr_w_t35(k,:)=k_cmr_w_t.';
%     
%     % dQ-Value3
%     md1_dq_w_t=fitlm(t(1:value3(4,k),1),dq_w_time(1:value3(4,k),1));
%     x_dq_w_t=md1_dq_w_t.Coefficients.Estimate;
%     b1lm_dq_w_t35(k,:)=x_dq_w_t.';
%     
%     md2_dq_w_t=fitlm(t(value3(4,k):60,1),dq_w_time(value3(4,k):60,1));
%     k_dq_w_t=md2_dq_w_t.Coefficients.Estimate;
%     b2lm_dq_w_t35(k,:)=k_dq_w_t.';
%     
%     % Q-Value4
%     md1_q_w_t=fitlm(t(1:value4(4,k),1),q_w_time(1:value4(4,k),1));
%     x_q_w_t=md1_q_w_t.Coefficients.Estimate;
%     b1lm_q_w_t35(k,:)=x_q_w_t.';
%     
%     md2_q_w_t=fitlm(t(value4(4,k):60,1),q_w_time(value4(4,k):60,1));
%     k_q_w_t=md2_q_w_t.Coefficients.Estimate;
%     b2lm_q_w_t35(k,:)=k_q_w_t.';
%     
%    
%     %find auc integrals up to plateau for different curves, i.e.,
%     %cbf, cmr, dqdt and q
%     
%     % Auc for cbf curve
%     auc_cbf35(k,1)=massi;
%     auc_cbf35(k,2)=trapz(t(1:value1(4,k),1),cbf_w_time(1:value1(4,k),1));
%     
%     % Auc for cmr curve
%     auc_cmr35(k,1)=massi;
%     auc_cmr35(k,2)=trapz(t(1:value2(4,k),1),cmr_w_time(1:value2(4,k),1));
%     
%     % Auc for dq curve
%     auc_dq35(k,1)=massi;
%     auc_dq35(k,2)=trapz(t(1:value3(4,k),1),cmr_w_time(1:value3(4,k),1));
%     
%     % Auc for q curve
%     auc_q35(k,1)=massi;
%     auc_q35(k,2)=trapz(t(1:value4(4,k),1),q_w_time(1:value4(4,k),1));
%  
%    

% 
% T_arterial = 36.3 %run through different arterial blood temperatures
%     [t,temph] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,37);
%     [t,tempi] = ode45(@(t,temp) eqn5_ischemic(t,temp,T_arterial),tspan,37);
%     cbf_w_timeh = zeros(length(t), 1);
%     cmr_w_timeh = zeros(length(t), 1);
%     dq_w_timeh = zeros(length(t), 1);
%     q_w_time = zeros(length(t),1); % added to calculate integral of dq/dt over t
%     cbf_w_timei = zeros(length(t), 1);
%     cmr_w_timei = zeros(length(t), 1);
%     dq_w_timei = zeros(length(t), 1);
% 
%     for i = 1:length(temph) % healthy tissue 
%         [cbf,cmr] = CBF_CMR_calculator(temph(i,1));
%         cbf_w_timeh(i,1) = cbf*(massh/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timeh(i,1) = cmr*10^6*(massh/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timeh(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temph(i,1)-T_arterial))*(massh/100); %kJ/min
%     end
%     
%     
%     for i = 1:length(tempi) % ischemic tissue
%         [normal_cbf, normal_cmr] = CBF_CMR_calculator(tempi(i,1));
%         % because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
%         % per 100 g if ASSUME density of brain = p_water = 1 g/ml)
%         if(normal_cbf >= 15)
%             cbf = 15;
%             cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
% %             temp=
%         else
%             cbf = normal_cbf;
%             %should w0 change? no
%             cmr = normal_cmr;
%         end
%         cbf_w_timei(i,1) = cbf*(massi/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timei(i,1) = cmr*10^6*(massi/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timei(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(tempi(i,1)-T_arterial))*(massi/100); %kJ/min
%     end
%     
%     cbf_w_time = cbf_w_timei + cbf_w_timeh;
%     cmr_w_time = cmr_w_timei + cmr_w_timeh;
%     dq_w_time = dq_w_timei + dq_w_timeh;
%     q_w_time = cumsum(dq_w_time);
%     
%     
%     figure(1)
%     plot(t,cbf_w_time);
%     hold on; 
%     
%     
%     figure(2)
%     plot(t,cmr_w_time);
%     hold on; 
%     
%     figure(3)
%     plot(t,dq_w_time);
%     hold on; 
%     
%     figure(4)
%     plot(t,q_w_time);
%     hold on
%     
%     value1(5,k)=find_plateau(cbf_w_time, t);
%     value2(5,k)=find_plateau(cmr_w_time, t);
%     value3(5,k)=find_plateau(dq_w_time, t);
%     value4(5,k)=find_plateau(q_w_time, t);
%     
%     % Rohit work-curve stripping
%     % CBF-Value1
%     md1_cbf_w_t=fitlm(t(1:value1(5,k),1),cbf_w_time(1:value1(5,k),1));
%     x_cbf_w_t=md1_cbf_w_t.Coefficients.Estimate;
%     b1lm_cbf_w_t36(k,:)=x_cbf_w_t.';
%     
%     md2_cbf_w_t=fitlm(t(value1(5,k):60,1),cbf_w_time(value1(5,k):60,1));
%     k_cbf_w_t=md2_cbf_w_t.Coefficients.Estimate;
%     b2lm_cbf_w_t36(k,:)=k_cbf_w_t.';
%     
%     % CMR-Value2
%     md1_cmr_w_t=fitlm(t(1:value2(5,k),1),cmr_w_time(1:value2(5,k),1));
%     x_cmr_w_t=md1_cmr_w_t.Coefficients.Estimate;
%     b1lm_cmr_w_t36(k,:)=x_cmr_w_t.';
%     
%     md2_cmr_w_t=fitlm(t(value2(5,k):60,1),cmr_w_time(value2(5,k):60,1));
%     k_cmr_w_t=md2_cmr_w_t.Coefficients.Estimate;
%     b2lm_cmr_w_t36(k,:)=k_cmr_w_t.';
%     
%     % dQ-Value3
%     md1_dq_w_t=fitlm(t(1:value3(5,k),1),dq_w_time(1:value3(5,k),1));
%     x_dq_w_t=md1_dq_w_t.Coefficients.Estimate;
%     b1lm_dq_w_t36(k,:)=x_dq_w_t.';
%     
%     md2_dq_w_t=fitlm(t(value3(5,k):60,1),dq_w_time(value3(5,k):60,1));
%     k_dq_w_t=md2_dq_w_t.Coefficients.Estimate;
%     b2lm_dq_w_t36(k,:)=k_dq_w_t.';
%     
%     % Q-Value4
%     md1_q_w_t=fitlm(t(1:value4(5,k),1),q_w_time(1:value4(5,k),1));
%     x_q_w_t=md1_q_w_t.Coefficients.Estimate;
%     b1lm_q_w_t36(k,:)=x_q_w_t.';
%     
%     md2_q_w_t=fitlm(t(value4(5,k):60,1),q_w_time(value4(5,k):60,1));
%     k_q_w_t=md2_q_w_t.Coefficients.Estimate;
%     b2lm_q_w_t36(k,:)=k_q_w_t.';
%     
%    
%     %find auc integrals up to plateau for different curves, i.e.,
%     %cbf, cmr, dqdt and q
%     
%     % Auc for cbf curve
%     auc_cbf36(k,1)=massi;
%     auc_cbf36(k,2)=trapz(t(1:value1(5,k),1),cbf_w_time(1:value1(5,k),1));
%     
%     % Auc for cmr curve
%     auc_cmr36(k,1)=massi;
%     auc_cmr36(k,2)=trapz(t(1:value2(5,k),1),cmr_w_time(1:value2(5,k),1));
%     
%     % Auc for dq curve
%     auc_dq36(k,1)=massi;
%     auc_dq36(k,2)=trapz(t(1:value3(5,k),1),cmr_w_time(1:value3(5,k),1));
%     
%     % Auc for q curve
%     auc_q36(k,1)=massi;
%     auc_q36(k,2)=trapz(t(1:value4(5,k),1),q_w_time(1:value4(5,k),1));
%  
%    
%  T_arterial = 37.3 %run through different arterial blood temperatures
%     [t,temph] = ode45(@(t,temp) eqn5(t,temp,T_arterial),tspan,37);
%     [t,tempi] = ode45(@(t,temp) eqn5_ischemic(t,temp,T_arterial),tspan,37);
%     cbf_w_timeh = zeros(length(t), 1);
%     cmr_w_timeh = zeros(length(t), 1);
%     dq_w_timeh = zeros(length(t), 1);
%     q_w_time = zeros(length(t),1); % added to calculate integral of dq/dt over t
%     cbf_w_timei = zeros(length(t), 1);
%     cmr_w_timei = zeros(length(t), 1);
%     dq_w_timei = zeros(length(t), 1);
% 
%     for i = 1:length(temph) % healthy tissue 
%         [cbf,cmr] = CBF_CMR_calculator(temph(i,1));
%         cbf_w_timeh(i,1) = cbf*(massh/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timeh(i,1) = cmr*10^6*(massh/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timeh(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temph(i,1)-T_arterial))*(massh/100); %kJ/min
%     end
%     
%     
%     for i = 1:length(tempi) % ischemic tissue
%         [normal_cbf, normal_cmr] = CBF_CMR_calculator(tempi(i,1));
%         % because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
%         % per 100 g if ASSUME density of brain = p_water = 1 g/ml)
%         if(normal_cbf >= 15)
%             cbf = 15;
%             cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
% %             temp=
%         else
%             cbf = normal_cbf;
%             %should w0 change? no
%             cmr = normal_cmr;
%         end
%         cbf_w_timei(i,1) = cbf*(massi/100); %convert cbf from ml/100g/min to ml/min
%         cmr_w_timei(i,1) = cmr*10^6*(massi/100); %convert cmr from mol02/100g/min to umol/min
%         dq_w_timei(i,1) = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(tempi(i,1)-T_arterial))*(massi/100); %kJ/min
%     end
%     
%     cbf_w_time = cbf_w_timei + cbf_w_timeh;
%     cmr_w_time = cmr_w_timei + cmr_w_timeh;
%     dq_w_time = dq_w_timei + dq_w_timeh;
%     q_w_time = cumsum(dq_w_time);
%     
%     
%     figure(1)
%     plot(t,cbf_w_time);
%     hold on; 
%     
%     
%     figure(2)
%     plot(t,cmr_w_time);
%     hold on; 
%     
%     figure(3)
%     plot(t,dq_w_time);
%     hold on; 
%     
%     figure(4)
%     plot(t,q_w_time);
%     hold on
%     
%     value1(6,k)=find_plateau(cbf_w_time, t);
%     value2(6,k)=find_plateau(cmr_w_time, t);
%     value3(6,k)=find_plateau(dq_w_time, t);
%     value4(6,k)=find_plateau(q_w_time, t);
%     
%     % Rohit work-curve stripping
%     % CBF-Value1
%     md1_cbf_w_t=fitlm(t(1:value1(6,k),1),cbf_w_time(1:value1(6,k),1));
%     x_cbf_w_t=md1_cbf_w_t.Coefficients.Estimate;
%     b1lm_cbf_w_t37(k,:)=x_cbf_w_t.';
%     
%     md2_cbf_w_t=fitlm(t(value1(6,k):60,1),cbf_w_time(value1(6,k):60,1));
%     k_cbf_w_t=md2_cbf_w_t.Coefficients.Estimate;
%     b2lm_cbf_w_t37(k,:)=k_cbf_w_t.';
%     
%     % CMR-Value2
%     md1_cmr_w_t=fitlm(t(1:value2(6,k),1),cmr_w_time(1:value2(6,k),1));
%     x_cmr_w_t=md1_cmr_w_t.Coefficients.Estimate;
%     b1lm_cmr_w_t37(k,:)=x_cmr_w_t.';
%     
%     md2_cmr_w_t=fitlm(t(value2(6,k):60,1),cmr_w_time(value2(6,k):60,1));
%     k_cmr_w_t=md2_cmr_w_t.Coefficients.Estimate;
%     b2lm_cmr_w_t37(k,:)=k_cmr_w_t.';
%     
%     % dQ-Value3
%     md1_dq_w_t=fitlm(t(1:value3(6,k),1),dq_w_time(1:value3(6,k),1));
%     x_dq_w_t=md1_dq_w_t.Coefficients.Estimate;
%     b1lm_dq_w_t37(k,:)=x_dq_w_t.';
%     
%     md2_dq_w_t=fitlm(t(value3(6,k):60,1),dq_w_time(value3(6,k):60,1));
%     k_dq_w_t=md2_dq_w_t.Coefficients.Estimate;
%     b2lm_dq_w_t37(k,:)=k_dq_w_t.';
%     
%     % Q-Value4
%     md1_q_w_t=fitlm(t(1:value4(6,k),1),q_w_time(1:value4(6,k),1));
%     x_q_w_t=md1_q_w_t.Coefficients.Estimate;
%     b1lm_q_w_t37(k,:)=x_q_w_t.';
%     
%     md2_q_w_t=fitlm(t(value4(6,k):60,1),q_w_time(value4(6,k):60,1));
%     
%     k_q_w_t=md2_q_w_t.Coefficients.Estimate;
%     b2lm_q_w_t37(k,:)=k_q_w_t.';
%     
%    
%     %find auc integrals up to plateau for different curves, i.e.,
%     %cbf, cmr, dqdt and q
%     
%     % Auc for cbf curve
%     auc_cbf37(k,1)=massi;
%     auc_cbf37(k,2)=trapz(t(1:value1(6,k),1),cbf_w_time(1:value1(6,k),1));
%     
%     % Auc for cmr curve
%     auc_cmr37(k,1)=massi;
%     auc_cmr37(k,2)=trapz(t(1:value2(6,k),1),cmr_w_time(1:value2(6,k),1));
%     
%     % Auc for dq curve
%     auc_dq37(k,1)=massi;
%     auc_dq37(k,2)=trapz(t(1:value3(6,k),1),cmr_w_time(1:value3(6,k),1));
%     
%     % Auc for q curve
%     auc_q37(k,1)=massi;
%     auc_q37(k,2)=trapz(t(1:value4(6,k),1),q_w_time(1:value4(6,k),1));
%  
   
 end

 
 
 



% figure(1)
% title('CBF vs Time for 300g Healthy Tissue + 200g Ischemic Tissue');
% ylabel('CBF (ml/min)');
% xlabel('time (min)');
% legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3');
% 
% figure(2)
% title('CMR vs Time for 300g Healthy Tissue + 200g Ischemic Tissue');
% ylabel('CMR (umol02/min)');
% xlabel('min');
% legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3');
% 
% figure(3)
% title('Heat Washout (dQ) vs Time for 300g Healthy Tissue + 200g Ischemic Tissue');
% ylabel('Rate of Change of Heat (kJ/min)');
% xlabel('min');
% legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3');
% figure(4)
% title('Heat Washout (Q) vs Time for 300g Healthy Tissue + 200g Ischemic Tissue');
% ylabel('Heat (kJ)');
% xlabel('min');
% legend('Perfusate temperature = 32.3','Perfusate temperature = 33.3', 'Perfusate temperature = 34.3', 'Perfusate temperature = 35.3', 'Perfusate temperature = 36.3', 'Perfusate temperature = 37.3');
% 
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
% disp('b1');
% disp(b1);
% disp('b2');
% disp(b2);



%%% Fitting lasso
% massh_pred=massh((:,1)