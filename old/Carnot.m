function [ P , carnot_eff ] = Carnot( T_h, T_l, flow, Cp)
%Carnot energy balance function 
% Assumptions: Infinte heat source & sink.
% Currently ingores heat exchangers, expanders, pumps, and condensors
% (everything)

% Inputs: Th -- high (source) temperature vector in Celsius 
%         Tl -- low (sink) temperature vector in Celsius
%         flow_h -- source flow rate vector in kg/s
%         Cp -- Specific heat (vector?) in J/(kg K) (4180 J/(kg K) for water)

% Outputs: 
%         P -- power (vector) extracted in W
%         carnot_eff -- fractional carnot efficiency

% convert C to K (when absolute T is needed, not relative)
Kelvin_h = T_h + 273;
Kelvin_l = T_l + 273;

% delta T -- difference in high and low temps (the same in C and K)
dT = T_h - T_l;

% Carnot efficiency 
carnot_eff = (1- Kelvin_l./Kelvin_h);

% Carnot power
P = Cp.*flow.*dT.*carnot_eff;

end

