function [myCBF, myCMR] = CBF_CMR_calculator_give_pars(temp, base_cbf, base_cmr)
%input temp in degrees celcius
%input base_cbf in ml/100g/min
%input base_cmr in molO2/100g/min
%do not use with temperatures below 25 deg C
%assumes region is half white matter, half gray matter(see w0)
alpha = 3; %unitless
beta = 0.1; %unitless
myCBF = base_cbf*alpha^(beta*(temp-37.3)); %returns CBF in ml/100g/min
myCMR = base_cmr*alpha^(beta*(temp-37.3)); %returns CMR02 in mol02/100g/min
