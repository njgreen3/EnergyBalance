function [ P_hy, T_or_H_o ] = IsentropicPower( T_or_H_i, p_i, p_o, m_dot, fluid, isen_eff, T_or_H )
% function [ P_hy ] = IsentropicPower( T_i, T_o, p_i, p_o, m_dot, fluid, isen_eff )
%EXPANDER Summary of this function
%   This process calculates the mechanical power provided or required by a 
%   generic expander device/system. Turbines and other expanders provide
%   mechanical power while pumps and compressors consume it. This is done
%   by calculating the change in enthalpy of fluid and multipying by the 
%   mass flow rate while entropy remains constant.  Non ideal processes are
%   taken into account through the isentropic efficiency.
% 
% Inputs:
%   T_or_H_i    Input temperature in Kelvin or enthalpy in J/kg
%   p_i         Input pressure of the fluid in Pa
%   p_o         Output temperature of the fluid in Pa
%   m_dot       Mass flow rate of the fluid in kg/s 
%   fluid       Name of the working fluid
%   isen_eff    Isentropic efficiency of the device. Expected range: 0-1
%               This input is optional. If it is not provided then it will 
%               be set to 1, which is the efficiency of an ideal isentropic
%               device where entropy is constant.
%               This value is treated differently depending on if there is
%               an increase or decrease in fluid enthalpy. If the enthalpy
%               decreases (energy is extracted from the system) then the
%               efficiency if multiplied by the change in energy because
%               less energy is extracted than the ideal case. If the
%               enthalpy increases (energy in added to the system) then the
%               efficiency is divided by the change in energy because more
%               energy is required than the ideal case.
%   T_or_H      String stating whether arguements 1&5 are temperature (T) 
%               or enthaply (H). Also applies to output
% 
% Outputs:
%   P_hy        Hydraulic power exctracted or consumed in W. Positive
%               values indicate the system is providing hydraulic power as
%               a turbine or expander. Negative values indicate the system 
%               is consuming hydraulic power as a compressor or pump.
%   T_o         Output temperature in Kelvin or enthalpy in J/kg

% This function makes use of the CoolProp wrapper for fluid properties. 
% The high level Matlab wrapper cannot handle vector inputs. There ways to
% make low level calls interfacing with the DLL so many values can be
% determined more quickly, but that is not implented here.
% More information about the Matlab wrapper can be found at
%   http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html
% A list of parameters can be found at 
%   http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table
% A list of fluids can be found at 
%   http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids

% % Check inputs
% if nargin < 7
% %   Set equal to 1 if isen_eff is not in input
%     isen_eff = 1;
%     
%     if nargin < 6
%         error(message('NotEnoughInputs'));
%     end
% end
% % check if isen_eff is in range
% if isen_eff > 1 || isen_eff < 0
%     error(message('isen_eff is out of range'));
% end

% Ensure temperature/enthalpy input variables are delt with appropriately
if upper(T_or_H) == 'T' %|| 'TEMPERATURE'
% Assume given temperature values are not at vaporization/condensation
    H_i = CoolProp.PropsSI('H','P',p_i,'T',T_or_H_i,fluid);
elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
% Enthalpy values can be in vaporization/condensation region    
    H_i = T_or_H_i;
else
    error('Please specify whether input uses variables of Temperature or Enthalpy.')
end

% Input enthalpy in J/kg and entropy in J/kg/K
% H_i = CoolProp.PropsSI('H', 'P', p_i, 'T', T_i, fluid);
S_i = CoolProp.PropsSI('S', 'P', p_i, 'H', H_i, fluid);

% Ideal Output entropy is equal to input entropy 
% Ideal Output enthalpy in J/kg and temperature
H_o_ideal = CoolProp.PropsSI('H', 'P', p_o, 'S', S_i, fluid);
% T_o_ideal = CoolProp.PropsSI('T', 'P', p_o, 'S', S_i, fluid);

% Calculate Ideal power consumed/provided
P_ideal = m_dot*(H_i - H_o_ideal);

% Check direction of power flow (into or out of system)
if P_ideal > 0
%   If positive then the system is providing power as a turbine
%   or expander.
%   This means the actual power provided is less than the ideal power
    P_hy = isen_eff*P_ideal;
elseif P_ideal < 0
%   If negative then the system is consuming power as a
%   compressor or pump.
%   This means the actual power consumed is greater than the ideal power.
    P_hy = P_ideal/isen_eff;
else
%   If the ideal power consumed/provided is 0 then either the fluid did not 
%   undergo any change or some (non-ideal) work was consumed to change it's
%   state
    warning('Actual mechanical power cannot be determined from the information provided.');
    P_hy = NaN;
    T_o = T_i;
    return;
end


H_o = H_i - P_hy/m_dot;
% T_o = CoolProp.PropsSI('T', 'P', p_o, 'H', H_o, fluid);

% Ensure temperature/enthalpy input variables are delt with appropriately
if upper(T_or_H) == 'T' %|| 'TEMPERATURE'
% Assume given temperature values are not at vaporization/condensation
    T_or_H_o = CoolProp.PropsSI('T', 'P', p_o, 'H', H_o, fluid);
elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
% Enthalpy values can be in vaporization/condensation region    
    T_or_H_o = H_o;
end
    
end

