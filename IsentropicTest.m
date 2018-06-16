clear variables
clc
addpath('C:\Users\njgreen3\Documents\MATLAB\CoolPropStuff')

% This test is for the IsentropicPower function. That function makes use of
% the high level CoolProp wrapper for Matlab. It can only hand individual
% input values and not vectors.

% The test will try compression and expansion scenarios.

%%
disp('This scenario compares an ideal and non-ideal expander under identical')
disp('conditions found in "Experimental testing of a low-temperature organic ')
disp('Rankine cycle (ORC) engine coupled with concentrating PV/thermal ')
disp('collectors: Laboratory and field tests" by Kosmadakis')
disp(' ')

Working_fluid = 'R404A'
T_or_H = 'T';
% Isentropic Efficiency-They did not report the expander's eff but I
% estimated a max based of other reported max values
is_eff = .8;
% Mass flow rate in kg/s
m_flow = 0.2237;

% temperature of input in degC
T_hot_in_C = 94.1;
% temperature of input in K
T_hot_in_K = T_hot_in_C + 273;

% temperature of output in degC
T_cool_out_C = 69.2;
% temperature of input in K
T_cool_out_K_reported = T_cool_out_C + 273;

% Pressure of input in bar
p_hot_in_bar = 28.6;
% Pressure of input in pa
p_hot_in_pa = p_hot_in_bar*100000;

% Pressure of input in bar
p_cool_out_bar = 13.9;
% Pressure of input in pa
p_cool_out_pa = p_cool_out_bar*100000;


disp('I calculated the following outputs in W and K:')
tic
[P_hy, T_cool_out_K] = IsentropicPower(T_hot_in_K, p_hot_in_pa, p_cool_out_pa, ...
    m_flow, Working_fluid, is_eff, T_or_H )
[P_ideal, T_cool_out_K_ideal] = IsentropicPower(T_hot_in_K, p_hot_in_pa, p_cool_out_pa, ...
    m_flow, Working_fluid, 1, T_or_H )
toc
tic
[P_hy, T_cool_out_K] = IsentropicPower_low(T_hot_in_K, p_hot_in_pa, p_cool_out_pa, ...
    m_flow, Working_fluid, is_eff, T_or_H )
[P_ideal, T_cool_out_K_ideal] = IsentropicPower_low(T_hot_in_K, p_hot_in_pa, p_cool_out_pa, ...
    m_flow, Working_fluid, 1, T_or_H )
toc
disp('The authors reported an electrical output of over 3 kW. This exceeds')
disp('my calcluated values when I would expect it to be smaller.')
disp('______________________________________________')
%%
disp('This scenario compares an ideal and non-ideal expander under identical')
disp('conditions found in "Design and experimental analysis of a mini ORC')
disp('power plant based on R245fa working fluid" by Galloni')
disp(' ')

Working_fluid = 'R245fa'
T_or_H = 'T';
% Isentropic Efficiency
is_eff = .849;
% Mass flow rate in kg/s
m_flow = 0.052;

% temperature of input to the expander in degC
T_hot_in_C = 89.7;
% temperature of input in K
T_hot_in_K = T_hot_in_C + 273;

% temperature of output to the expander in degC
T_cool_out_C = 38.4;
% temperature of input in K
T_cool_out_K_reported = T_cool_out_C + 273

% Pressure of input in bar
p_hot_in_bar = 9.95;
% Pressure of input in pa
p_hot_in_pa = p_hot_in_bar*100000;

% Pressure of input in bar
p_cool_out_bar = 2.02;
% Pressure of input in pa
p_cool_out_pa = p_cool_out_bar*100000;

disp('I calculated the following outputs in W and K:')
[P_hy, T_cool_out_K] = IsentropicPower(T_hot_in_K, p_hot_in_pa, p_cool_out_pa, ...
    m_flow, Working_fluid, is_eff, T_or_H )
[P_ideal, T_cool_out_K_ideal] = IsentropicPower(T_hot_in_K, p_hot_in_pa, p_cool_out_pa, ...
    m_flow, Working_fluid, 1, T_or_H )
disp('The authors reported an adiabatic output of 1.17 kW. This is less than')
disp('my calcluated values. I would expect it to be comparable to my ideal value.')
disp('______________________________________________')

%%
disp('This scenario compares an ideal and non-ideal pump under identical')
disp('conditions found in "Operation and performance of a low temperature')
disp('organic Rankine cycle" by Miao')
disp(' ')

Working_fluid = 'R123'
T_or_H = 'T';
% Isentropic Efficiency
is_eff = .74;
% Mass flow rate in kg/s
m_flow = 0.081;

% temperature of input to the expander in degC
T_in_C = 18.9;
% temperature of input in K
T_in_K = T_in_C + 273;

% temperature of output to the expander in degC
T_out_C = 18.9;
% temperature of input in K
T_out_K_reported = T_out_C + 273

% Pressure of input in bar
p_in_bar = 1.35;
% Pressure of input in pa
p_in_pa = p_in_bar*100000;

% Pressure of input in bar
p_out_bar = 10.3;
% Pressure of input in pa
p_out_pa = p_out_bar*100000;

disp('I calculated the following outputs in W and K:')
[P_hy, T_out_K] = IsentropicPower(T_in_K, p_in_pa, p_out_pa, ...
    m_flow, Working_fluid, is_eff,T_or_H )
[P_ideal, T_out_K_ideal] = IsentropicPower(T_in_K, p_in_pa, p_out_pa, ...
    m_flow, Working_fluid, 1, T_or_H )
disp('The authors reported an adiabatic consumption of 0.45 kW. This is much more than')
disp('my calcluated values. I would expect it to be comparable to my ideal value.')
disp('The discrepancy is probably caused by my choise of pressure values. The')
disp('paper only reports the maximum and minium pressure experienced by the')
disp('working fluid. I assumed the only singificant pressure changes occured at')
disp('the pump and the expander, but that is probably not the case. I also assumed')
disp('no temperature change at the pump, but I expect the pressure assumption')
disp('had a bigger impact on the discrepency.')
disp('____________________________________________________')