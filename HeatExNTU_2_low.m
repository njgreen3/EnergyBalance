function [ P, Tho, Tco ] = HeatExNTU_2_low( Thi, m_dot_h, p_h, Tci, m_dot_c, p_c, fluid_h, ...
fluid_c, A, U, HeatExType)
%HeatExNTU Summary 
%   This function calculates flow and outlet temperatures of a heat exchanger
%   using the Number of Transfer Units (NTU) method as described in
%   Fundementals of Heat and Mass Transfer by F. P. Incropera, D. P.
%   DeWitt, T. L. Bergman & A. S. Lavine. The method can be used when only
%   the inlet tempertures are known as opposed to inlet and outlet
%   temperatures.
%
%   This is done by calculating effectiveness which is defined as the ratio
%   of actual heat transferred to the maximum possible heat transferred for 
%   an infinite length counter-flow heat exchanger with the given inlet
%   temperatures. The effectivenes is calculated differently depending on 
%   the geometry of exchanger but is always a function of the ratio of heat
%   capacity rates as well as the number of transfer units (NTU). The NTU
%   depends on the overall heat transfer coefficient, the heat transfer
%   area, and the smaller of the two heat capacity rates. 
%   
%   In the event of a phase change the rate of change of energy required to
%   get the fluid just beyond the point of vaporizaion or condensation is
%   calculated using the CoolProp function. Then the function is called
%   again at the new temperature and the final temperatures and heat rates
%   are calucluated using a reduced heat exhange area. I currently assume
%   the effected heat transfer area is a half of the original input. This 
%   function cannot accurately handle both fluids changing phase because it
%   does not check which fluid would change first.
% 
% Note: The effectiveness is separte from efficiency. The calculations
%   assumes an efficiency of 1 where all heat is transferred from the hot
%   fluid to the cool fluid or remains within the the fluid. This also
%   assumes no fluid is undergoing a change in state, which may not be
%   sufficient for my final model.
%
% This function makes use of CoolProp, a shared DLL. 
% More information about the Matlab wrapper can be found at
%   http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html
% A list of parameters can be found at 
%   http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table
% A list of fluids can be found at 
%   http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
%
% Inputs:
%   Thi         Hot fluid inlet temperature in Kelvin
%   m_dot_h     Hot fluid mass flow rate in kg/s
%   p_h         Hot fluid pressure in pa
%   fluid_h     Name of Hot fluid (string)
%   Tci         Cool fluid inlet temperature in Kelvin
%   m_dot_c     Cool fluid mass flow rate in kg/s
%   p_c         Cool fluid pressure in pa
%   fluid_c     Name of Cool fluid (string)
%   HeatExType  String describing the heat exchanger geometry. Currently
%               this variable only accepts 'counter' and 'parallel', but
%               can be expanded to other exchanger/flow types
%               'counter' for counter flow
%               'parallel' for parallel
%   U           Overall heat transfer coefficient in W/(K m2)
%   A           Heat transfer area in m2
%
% Output:
%   P           Rate of heat exchenged (power) in W
%   Tho         Hot fluid outlet temperature in input temp units (K)
%   Tco         Cool fluid outlet temperature in input tempt units (K)

% initailize certain parameters used to call coolprop library
buffer_size = 1000;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);
backend = 'HEOS';
% coder.extrinsic('calllib') 
% coder.extrinsic('libpointer') 

% All the input signals should be the same size, so it shouldn't matter
% which one is used to obtain signal length
var_length = length(Thi);

% initialize handles for hot and cool fluids
[handle_h, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_h,ierr,herr,buffer_size);
[handle_c, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_c,ierr,herr,buffer_size);

% determine input parameter pair indecies for coolprop calls
[pq_ind,~] = calllib('coolprop','get_input_pair_index','PQ_INPUTS'); %pressure and Mass vapor quality
[pt_ind,~] = calllib('coolprop','get_input_pair_index','PT_INPUTS'); %pressure and temperature
[hp_ind,~] = calllib('coolprop','get_input_pair_index','HmassP_INPUTS'); %Mass specific enthalpy and pressure

% determine output parameter indecies for coolprop calls
[t_ind, ~] = calllib('coolprop','get_param_index','T'); %temperature
[h_ind, ~] = calllib('coolprop','get_param_index','Hmass'); %Mass specific enthalpy
[c_ind, ~] = calllib('coolprop','get_param_index','Cpmass'); %Mass specific constant pressure specific heat

% find pointer location of pressure and temperature input values
p_h_Ptr = libpointer('doublePtr',p_h);
p_c_Ptr = libpointer('doublePtr',p_c);
Thi_Ptr = libpointer('doublePtr',Thi);
Tci_Ptr = libpointer('doublePtr',Tci);

% Determine vaporization temperature for fluids at respective pressures.
% Initialize pointer locations for relevant variables
Q_Ptr = libpointer('doublePtr', zeros(var_length,1)); % Q is specifically set to 0
Tvh = zeros(var_length,1); % initialized to 0 temporarily
Tvc = zeros(var_length,1); % initialized to 0 temporarily
Tvh_Ptr = libpointer('doublePtr', Tvh);
Tvc_Ptr = libpointer('doublePtr', Tvc);

% Update values located at the Tvc & Tvh pointer locations for vaporization
% conditions
calllib('coolprop','AbstractState_update_and_1_out', handle_h, pq_ind, p_h_Ptr, Q_Ptr, var_length, t_ind, Tvh_Ptr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', handle_c, pq_ind, p_c_Ptr, Q_Ptr, var_length, t_ind, Tvc_Ptr, ierr,herr,buffer_size);

% Set Tvh & Tvc to values stored at pointer location
Tvh = get(Tvh_Ptr,'Value');
Tvc = get(Tvc_Ptr,'Value');

% Determine the specific heats of the incoming fluids using coolprop
% Initialize pointer locations for relevant variables
Ch = zeros(var_length,1); % initialized to 0 temporarily
Cc = zeros(var_length,1); % initialized to 0 temporarily
Ch_Ptr = libpointer('doublePtr', Ch); 
Cc_Ptr = libpointer('doublePtr', Cc); 

% Update values located at the Ch & Cc pointer locations for input
% conditions
calllib('coolprop','AbstractState_update_and_1_out', handle_h, pt_ind, p_h_Ptr, Thi_Ptr, var_length, c_ind, Ch_Ptr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', handle_c, pt_ind, p_c_Ptr, Tci_Ptr, var_length, c_ind, Cc_Ptr, ierr,herr,buffer_size);

% Set Ch & Cc to values stored at pointer location
Ch = get(Ch_Ptr,'Value');   % J/(kg K)
Cc = get(Cc_Ptr,'Value');   % J/(kg K)

% If the fluid is withing 1e-4 units of the vaporization point then the 
% will return an inaccurate values. If that happens then the specific heat 
% will be set to inf so its temperature will not change even though heat is
% exchanged.
Ch(abs(Thi - Tvh) < 1e-4) = inf;    % J/(kg K)
Cc(abs(Tci - Tvc) < 1e-4) = inf;    % J/(kg K)

% Calculate the heat capacity rates of the hot and cool fluids
Ch_rate = Ch .* m_dot_h;     % W/K
Cc_rate = Cc .* m_dot_c;     % W/K

% Determine the maximum and minimum of the two heat capacity rates
Cmax = Ch_rate; % W/K
Cmax(Ch_rate < Cc_rate) = Cc_rate(Ch_rate < Cc_rate);   % W/K
Cmin = Ch_rate; % W/K
Cmin(Ch_rate > Cc_rate) = Cc_rate(Ch_rate > Cc_rate);   % W/K


% Calcluate the heat capacity rate ratio 
Cr = Cmin./Cmax;

% Calculate the number of transfer units of the heat exchanger
NTU = U * A ./ Cmin;

% Calculate effectiveness based on type of heat exchanger
switch HeatExType
    case 'counter'  % If the heat exchanger is counter-flow
        effectiveness = (1 - exp(-NTU.*(1 - Cr))) ./ ...
                        (1 - Cr.*exp(-NTU.*(1 - Cr)));
    case 'parallel' % If the heat exchanger is parallel-flow
        effectiveness = (1 - exp(-NTU.*(1 + Cr))) ./ ...
                        (1 + Cr);
    otherwise       % If not exchanger specified
%         warning('No heat exchanger type specified, therefore no heat is exchanged.')
        effectiveness = 0;
end        

% % However if the Cr is 0 (Cmax is inf) then effectiveness is calcluated
% % regardless of the heat exchanger flow.
% effectiveness(Cr == 0) =  1 - exp(-NTU(Cr == 0));

% Maximum energy transfer rate is equal for an infinite length heat
% exchanger is equal to Cmin multiplied with the difference in inlet
% temperatures. The actual energy transfer rate is the maximum rate
% multiplied by the effectiveness (Assuming no heat loss).
P = effectiveness.*Cmin.*(Thi - Tci);

% The change in temperature is calculated for both sides by dividing the
% energy tranfer rate by the heat capacity rates for both sides. This
% assumes the energy transfer rate is the same for both sides.
deltaTh = P./Ch_rate;
deltaTc = P./Cc_rate;

% The outlet temperatures are caluctaed by subtracting the change in
% temperature from the hot fluid inlet temperature and adding the change in
% temperature from the cool fluid inlet temperature.
Tho = Thi - deltaTh;
Tco = Tci + deltaTc;
    
% check if hot fluid would condense 
if any(Tho < Tvh & Thi > Tvh)
%    disp('Condensation')
%  Determine the enthalpy extracted from the hot fluid to condense it from 
%  starting temeprature. A slight 1e-4 K temperature change beyond
%  condensation is subtracted because coolprop cannot handle fluids within 
%  1e-4 K of the vaporization/condensation curve.
   Hh_init = zeros(var_length,1); % initialized to 0 temporarily
   Hh_cond = zeros(var_length,1); % initialized to 0 temporarily
   Hh_init_Ptr = libpointer('doublePtr', Hh_init); 
   Hh_cond_Ptr = libpointer('doublePtr', Hh_cond);
   
%  Determine value and point for a temperature just shy of the hot fluid
%  vaporiztion temperature
   Tvh_min = Tvh - 1e-4;
   Tvh_min_Ptr = libpointer('doublePtr', Tvh_min); 

% Update values located at the Hh_init & Hh_cond pointer locations for input
% conditions
   calllib('coolprop','AbstractState_update_and_1_out', handle_h, pt_ind, p_h_Ptr, Thi_Ptr, var_length, h_ind, Hh_init_Ptr, ierr,herr,buffer_size);
   calllib('coolprop','AbstractState_update_and_1_out', handle_h, pt_ind, p_h_Ptr, Tvh_min_Ptr, var_length, h_ind, Hh_cond_Ptr, ierr,herr,buffer_size);
   
% Set H to values stored at pointer locations
   Hh_init = get(Hh_init_Ptr,'Value');   
   Hh_cond = get(Hh_cond_Ptr,'Value');   

%  Calculate rate of heat transfer due to change in temp and condensation
   P_cond = (Hh_init - Hh_cond).*m_dot_h;
   
%  recalucluate the change in temp for the cool fluid due to condensation
%  of the hot fluid
   deltaTc = P_cond./Cc_rate;
%  recalculate the temperature of the cool fluid 
   Tc_cond = Tci + deltaTc;
   
%  Run this function again but with the input temps set to temps after
%  condensation. Because the caluclation for the heat exchanger is broken
%  up to multiple parts, the effective area should be reduced. I will
%  temporarily assume it is reduced by a factor of two but should return to
%  this later.
   [ P_deltaT, Tho, Tco ] = HeatExNTU_2_low( Tvh_min, m_dot_h, p_h, ...
       Tc_cond, m_dot_c, p_c, fluid_h, fluid_c, A/2, U, HeatExType);
   P = P_cond + P_deltaT;

% else check if cool fluid would vaporize
elseif any(Tco > Tvc & Tci < Tvc)
%    disp('Vaporization')
%  Determine the enthalpy needed to vaporize the cool fluid from 
%  starting temeprature. A slight 1e-4 K temperature change beyond
%  vaporization is added because coolprop cannot handle fluids within 
%  1e-4 K of the vaporization/condensation curve.
   Hc_init = zeros(var_length,1); % initialized to 0 temporarily
   Hc_vap = zeros(var_length,1); % initialized to 0 temporarily
   Hc_init_Ptr = libpointer('doublePtr', Hc_init); 
   Hc_vap_Ptr = libpointer('doublePtr', Hc_vap); 

%  Determine value and point for a temperature just over the cool fluid
%  vaporiztion temperature
   Tvc_max = Tvh + 1e-4;
   Tvc_max_Ptr = libpointer('doublePtr', Tvc_max); 
   
% Update values located at the Hc_init & Hc_vap pointer locations for input
% conditions
   calllib('coolprop','AbstractState_update_and_1_out', handle_c, pt_ind, p_c_Ptr, Tci_Ptr, var_length, h_ind, Hc_init_Ptr, ierr,herr,buffer_size);
   calllib('coolprop','AbstractState_update_and_1_out', handle_c, pt_ind, p_c_Ptr, Tvc_max_Ptr, var_length, h_ind, Hc_vap_Ptr, ierr,herr,buffer_size);

% Set H to values stored at pointer locations
   Hc_init = get(Hc_init_Ptr,'Value');   
   Hc_vap = get(Hc_vap_Ptr,'Value');   
     
%  Calculate rate of heat transfer due to change in temp and vaporization
   P_vap = (Hc_vap - Hc_init).*m_dot_c;
   
%  recalucluate the change in temp for the hot fluid due to vaporzation of
%  the cool fluid
   deltaTh = P_vap./Ch_rate;
%  recalculate the temperature of the hot fluid 
   Th_vap = Thi - deltaTh;
   
%  Run this function again but with the input temps set to temps after
%  vaporization. Because the caluclation for the heat exchanger is broken
%  up to multiple parts, the effective area should be reduced. I will
%  temporarily assume it is reduced by a factor of two but should return to
%  this later.
   [ P_deltaT, Tho, Tco ] = HeatExNTU_2_low( Th_vap, m_dot_h, p_h, ...
       Tvc + 1e-4, m_dot_c, p_c, fluid_h, fluid_c, A/2, U, HeatExType);
   P = P_vap + P_deltaT;
end

% check if hot fluid would condense and cool fluid would vaporize
if any((Tho < Tvh & Thi > Tvh) & (Tco > Tvc & Tci < Tvc))
    error('This function cannot accurately handle both fluids changing phase.')
end

end

