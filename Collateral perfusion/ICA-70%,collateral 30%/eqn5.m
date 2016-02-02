function dtemp = eqn5(t, temp, T_arterial)
% computes the RHS of the bioheat equation as formulated in Yablonskiy's
% paper. When passed an arterial blood temperature and solved, will return
% a function temp(t) where t is time.
% The dynamic bio-heat equation describes the rate of temperature change
% in the brain, temp', when the resting state is disturbed by global changes
% in blood flow, incoming blood temperature, or oxygen consumption
H0= 470; %  kJ/mol O2; from Yablonskiy paper,2000
Hb= 28; % kJ/mol O2; from Yablonskiy paper
[cbf, cmr] = CBF_CMR_calculator(temp);
p_blood = 1; % g/ml; assumed to be same as for water, Yablonskiy
c_blood = 4.178*10^-3; % specific heat in kJ/g/(degree celcius change); assumed to be same as for water, Yablonskiy
% brain mass? brain volume = 500 ml, assuming p_brain = 1 g/ml, mass = 500 g
c_tissue = .37; % specific heat, kJ/(delC*100g) taken from Konstas table 4, equal for white and gray

dtemp = 0.7*((H0-Hb)*cmr - p_blood*c_blood*cbf*(temp-T_arterial))/c_tissue; %dtemp/dt = delC/min