function [myCBF, myCMR] = CBF_CMR_calculator(temp)
%input temp in degrees celcius
%do not use with temperatures below 25 deg C
%assumes region is half white matter, half gray matter(see w0)
alpha = 3; %unitless
beta = 0.1; %unitless
w0 = 90; %ml/100g/min; average between gray and white matter; Konstas paper
q0 = 220*10^(-6); %mol O2/100 g brain tissue/min; taken from Yablonskiy paper, but changed units
myCBF = w0*alpha^(beta*(temp-38.3)); %returns CBF in ml/100g/min
myCMR = q0*alpha^(beta*(temp-38.3)); %returns CMR02 in mol02/100g/min
