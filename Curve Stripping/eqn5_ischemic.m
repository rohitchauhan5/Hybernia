function dtemp = eqn5_ischemic(t, temp, T_arterial)
% computes the RHS of the bioheat equation as formulated in Yablonskiy's
% paper. When passed an arterial blood temperature and solved, will return
% a function temp(t) where t is time.
% The dynamic bio-heat equation describes the rate of temperature change
% in the brain, temp', when the resting state is disturbed by global changes
% in blood flow, incoming blood temperature, or oxygen consumption
H0= 470; %  kJ/mol O2; from Yablonskiy paper,2000
Hb= 28; % kJ/mol O2; from Yablonskiy paper
[normal_cbf, normal_cmr] = CBF_CMR_calculator(temp);
% because ischemic, assuming CBF = 25 ml blood/100 g tissue/min (can say
% per 100 g if ASSUME density of brain = p_water = 1 g/ml)
    if(normal_cbf >= 15)
        cbf = 15;
        cmr = cbf/normal_cbf*normal_cmr; %1:1 change in cbf & cmr
    else
        cbf = normal_cbf;
        cmr = normal_cmr;
    end
    
p_blood = 1; % g/ml; assumed to be same as for water, Yablonskiy
c_blood = 4.178*10^-3; % specific heat in kJ/g/(degree celcius change); assumed to be same as for water, Yablonskiy
c_tissue = .37; % kJ/100 g/(deg Celcius change); taken from Konstas table 4; equal for white and gray
%((H0-Hb)*cmr - p_blood*c_blood*cbf*(temp-T_arterial))/c_tissue
dtemp = ((H0-Hb)*cmr - p_blood*c_blood*cbf*(temp-T_arterial))/c_tissue;
