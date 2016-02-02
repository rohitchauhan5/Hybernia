cbf_vector = zeros(1,length(24:38));
cmr_vector = zeros(1,length(24:38));
temp = 24:38;
for j = 0:10
    for temp1 = 24:38 
        [normal_cbf, normal_cmr] = CBF_CMR_calculator(temp1);
        % because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
        % per 100 g if ASSUME density of brain = p_water = 1 g/ml)
        if(normal_cbf >= 35)
        cbf = 35*0.1*j+normal_cbf*0.1*(10-j);
        cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
        else
        cbf = normal_cbf;
        %should w0 change? no
        cmr = normal_cmr;
        end
        cbf_vector(1,temp1-23) = cbf;
        cmr_vector(1,temp1-23) = cmr*10^6; % CMR in umol02/100g/min
    
    end
    figure(3)
    plot(temp, cbf_vector);
    hold on;
    figure(4)
    plot(temp, cmr_vector);
    hold on;
end
%% plot cbf vs temp
figure(3)
plot(temp, cbf_vector);
axis([24, 38, 0, 60]);
title('CBF in ischemic tissue');
xlabel('temperature (deg Celcius)');
ylabel('CBF (ml/100g/min)');
set(gca,'XDir','Reverse')
%print('plot_cbf_ischemic', '-djpeg');

%% plot cmr vs temp
figure(4)
plot(temp, cmr_vector);
axis([24, 38, 30,180]);
title('CMR in ischemic tissue');
xlabel('temperature (deg Celcius)');
ylabel('CMR (umol02/100g/min)');
set(gca,'XDir','Reverse')
%print('plot_cmr_ischemic', '-djpeg');
