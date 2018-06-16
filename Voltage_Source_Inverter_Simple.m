function [ Pdc, Ploss ] = Voltage_Source_Inverter_Simple(Pload, eff)
% function [ Pdc, Ploss ] = Voltage_Source_Inverter_Simple(Pload, Qload, Vdc, Vac, eff, method_flag )
                            
%INVERTER Summary of this function goes here
%   This function simulates a very simple three phase voltage source 
%   inverter for an energy balance microgrid model. For a given load, input
%   and output voltage levels, and efficiency, the scripts will calculated
%   the power comsumptions seen by the DC source as well as the power loss
%   in the inverter.

if nargin == 0
   Pload = 1000;    % W
%    Qload = 400;    % VAR
%    Vdc = 480;   % V
%    Vac = 480;   % V
   eff = 0.95;   % unitless
%    method_flag = 1;
end

% if method_flag
% first method
% Active Power out = Active Power in * efficiency
    Pdc = Pload./eff;

% else
% % Alternate method
% % |Apparent Power out| = |Apparent Power in * efficiency|
%     Sload = abs(Pload + j*Qload);
%     Pdc = Sload/eff;
% end

% loss is the difference in power
Ploss = Pdc - Pload;

end

