function [my_dq] = dq_calculator(temp)
%input temp in degrees celcius
%do not use with temperatures below 25 deg C
%assumes region is half white matter, half gray matter(see w0)
% alpha = 3; %unitless
% beta = 0.1; %unitless
% w0 = 50; %ml/100g/min; average between gray and white matter; Konstas paper
% q0 = 1.5*10^(-4); %mol O2/100 g brain tissue/min; taken from Yablonskiy paper, but changed units
c_tissue = .37; % specific heat, kJ/(delC*100g) taken from Konstas table 4, equal for white and gray
my_dq= c_tissue*temp; %returns CBF in ml/100g/min
% myCMR = q0*alpha^(beta*(temp-37)); %returns CMR02 in mol02/100g/min
