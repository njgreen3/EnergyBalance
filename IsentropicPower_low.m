function [ P_hy, T_or_H_o ] = IsentropicPower_low( T_or_H_i, p_i, p_o, m_dot, fluid, isen_eff, T_or_H )
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

% initailize certain parameters used to call coolprop library
buffer_size = 255;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);

coder.extrinsic('calllib') 
coder.extrinsic('libpointer') 

% All the input signals should be the same size, so it shouldn't matter
% which one is used to obtain signal length
var_length = length(T_or_H_i);
% T_or_H_i

fluid = 'R245fa';
% initialize handles for hot and cool fluids
[handle, ~] = calllib('coolprop','AbstractState_factory',backend,fluid,ierr,herr,buffer_size);

% determine input parameter pair indecies for coolprop calls
% [pq_ind,~] = calllib('coolprop','get_input_pair_index','PQ_INPUTS'); %pressure and Mass vapor quality
[pt_ind,~] = calllib('coolprop','get_input_pair_index','PT_INPUTS'); %pressure and temperature
[hp_ind,~] = calllib('coolprop','get_input_pair_index','HmassP_INPUTS'); %Mass specific enthalpy and pressure
[ps_ind,~] = calllib('coolprop','get_input_pair_index','PSmass_INPUTS'); %pressure and Mass specific entropy

% determine output parameter indecies for coolprop calls
[t_ind, ~] = calllib('coolprop','get_param_index','T'); %temperature
[s_ind, ~] = calllib('coolprop','get_param_index','S'); %Mass specific entropy
[h_ind, ~] = calllib('coolprop','get_param_index','Hmass'); %Mass specific enthalpy
% [c_ind, ~] = calllib('coolprop','get_param_index','Cpmass'); %Mass specific constant pressure specific heat

% find pointer location of pressure input and output values
p_i_Ptr = libpointer('doublePtr',p_i);
p_o_Ptr = libpointer('doublePtr',p_o);

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
% Assume given temperature values are not at vaporization/condensation.
% Initialize T and H inputs
    T_i = T_or_H_i;
    T_i_Ptr = libpointer('doublePtr',T_i);
    H_i = zeros(var_length,1); % initialized to 0 temporarily
    H_i_Ptr = libpointer('doublePtr',H_i);
% Update values located at Hhi and Hci pointer locations
    calllib('coolprop','AbstractState_update_and_1_out', handle, pt_ind, p_i_Ptr, T_i_Ptr, var_length, h_ind, H_i_Ptr, ierr,herr,buffer_size);
% Set Hhi & Hci to values stored at pointer location
    H_i = get(H_i_Ptr,'Value');
    
elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
% Enthalpy values can be in vaporization/condensation region  
% Initialize H values
    H_i = T_or_H_i;
    H_i_Ptr = libpointer('doublePtr',H_i);
else
    error('Please specify whether input uses variables of Temperature or Enthalpy.')
end

% Input entropy in J/kg/K
% Initialize pointer locations for relevant variables
S_i_Ptr = libpointer('doublePtr', zeros(var_length,1)); 
% Update values at S_i pointer location
calllib('coolprop','AbstractState_update_and_1_out', handle, hp_ind, H_i_Ptr, p_i_Ptr, var_length, s_ind, S_i_Ptr, ierr,herr,buffer_size);
% Set S_i to values stored at pointer location
% S_i = get(S_i_Ptr,'Value');

% Ideal Output entropy is equal to input entropy 
% Ideal Output enthalpy in J/kg and temperature
H_o_ideal_Ptr = libpointer('doublePtr',zeros(var_length,1));
% Update values located at Hhi and Hci pointer locations
calllib('coolprop','AbstractState_update_and_1_out', handle, ps_ind, p_o_Ptr, S_i_Ptr, var_length, h_ind, H_o_ideal_Ptr, ierr,herr,buffer_size);
% Set H_o_ideal to values stored at pointer location
H_o_ideal = get(H_o_ideal_Ptr,'Value');
% T_o_ideal = CoolProp.PropsSI('T', 'P', p_o, 'S', S_i, fluid);

% Calculate Ideal power consumed/provided
P_ideal = m_dot.*(H_i - H_o_ideal);

% Initialize hydrolic power variable
P_hy = zeros(var_length,1);
%   If positive then the system is providing power as a turbine
%   or expander.
%   This means the actual power provided is less than the ideal power
P_hy(P_ideal > 0) = isen_eff.*P_ideal(P_ideal > 0);

%   If negative then the system is consuming power as a
%   compressor or pump.
%   This means the actual power consumed is greater than the ideal power.
P_hy(P_ideal < 0) = P_ideal(P_ideal < 0)./isen_eff;

% %   If the ideal power consumed/provided is 0 then either the fluid did not 
% %   undergo any change or some (non-ideal) work was consumed to change it's
% %   state
%     warning('Actual mechanical power cannot be determined from the information provided.');
%     P_hy = NaN;
%     T_o = T_i;
%     return;


H_o = H_i - P_hy./m_dot;
H_o_Ptr = libpointer('doublePtr',H_o);
% T_o = CoolProp.PropsSI('T', 'P', p_o, 'H', H_o, fluid);

% Ensure temperature/enthalpy input variables are delt with appropriately
if upper(T_or_H) == 'T' %|| 'TEMPERATURE'
% Assume given temperature values are not at vaporization/condensation
    T_or_H_o_Ptr = libpointer('doublePtr',zeros(var_length,1));
    calllib('coolprop','AbstractState_update_and_1_out', handle, hp_ind, H_o_Ptr, p_o_Ptr, var_length, t_ind, T_or_H_o_Ptr, ierr,herr,buffer_size);
    T_or_H_o = get(T_or_H_o_Ptr,'Value');
    
elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
% Enthalpy values can be in vaporization/condensation region    
    T_or_H_o = H_o;
end
calllib('coolprop','AbstractState_free', handle, ierr,herr,buffer_size);
end

