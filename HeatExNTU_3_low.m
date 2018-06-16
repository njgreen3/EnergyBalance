function [ P, T_or_H_ho, T_or_H_co ] = HeatExNTU_3_low( T_or_H_hi, m_dot_h, p_h,  ...
    T_or_H_ci, m_dot_c, p_c, fluid_h, fluid_c, HeatExType, U, A, T_or_H )
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
%   T_or_H_hi   Hot fluid inlet temperature in Kelvin or enthalpy in J/kg
%   m_dot_h     Hot fluid mass flow rate in kg/s
%   p_h         Hot fluid pressure in pa
%   fluid_h     Name of Hot fluid (string)
%   T_or_H_ci   Cool fluid inlet temperature in Kelvin or enthalpy in J/kg
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
%   T_or_H      String stating whether arguements 1&5 are temperature (T) 
%               or enthaply (H). Also applies to output
%
% Output:
%   P           Rate of heat exchenged (power) in W
%   T_or_H_ho   Hot fluid outlet temperature (K) or enthalpy J/kg
%   T_or_H_co   Cool fluid outlet temperature (K) or enthalpy J/kg

% initailize certain parameters used to call coolprop library
buffer_size = 255;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);
coder.extrinsic('calllib') 
coder.extrinsic('libpointer') 

% All the input signals should be the same size, so it shouldn't matter
% which one is used to obtain signal length
var_length = length(T_or_H_hi);
% T_or_H_hi

% fluid_h = 'water';
% fluid_c = 'R245fa';
% initialize handles for hot and cool fluids
% [fluid_h, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_h,ierr,herr,buffer_size);
% [fluid_c, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_c,ierr,herr,buffer_size);

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
% Thi_Ptr
% Tci_Ptr
% Hhi_Ptr
% Hci_Ptr

% initialize output to zeros of appropriate length. It is necessary to
% assign output variables before recursively calling this function.
T_or_H_ho = zeros(var_length,1);
T_or_H_co = zeros(var_length,1);
% P = zeros(size(T_or_H_hi));

% T_or_H = 'T';
% Ensure temperature/enthalpy input variables are delt with appropriately
% if upper(T_or_H) == 'T' %|| 'TEMPERATURE'
if T_or_H == 'T' %|| 'TEMPERATURE'
% Assume given temperature values are not at vaporization/condensation

% Initialize T and H values 
    Thi = T_or_H_hi;
    Tci = T_or_H_ci;
    Thi_Ptr = libpointer('doublePtr',Thi);
    Tci_Ptr = libpointer('doublePtr',Tci);
    Hhi = zeros(var_length,1); % initialized to 0 temporarily
    Hci = zeros(var_length,1); % initialized to 0 temporarily
    Hhi_Ptr = libpointer('doublePtr',Hhi);
    Hci_Ptr = libpointer('doublePtr',Hci);
    
% Update values located at Hhi and Hci pointer locations
    calllib('coolprop','AbstractState_update_and_1_out', fluid_h, pt_ind, p_h_Ptr, Thi_Ptr, var_length, h_ind, Hhi_Ptr, ierr,herr,buffer_size);
    calllib('coolprop','AbstractState_update_and_1_out', fluid_c, pt_ind, p_c_Ptr, Tci_Ptr, var_length, h_ind, Hci_Ptr, ierr,herr,buffer_size);
% Set Hhi & Hci to values stored at pointer location
    Hhi = get(Hhi_Ptr,'Value');
    Hci = get(Hci_Ptr,'Value');

% elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
elseif T_or_H == 'H' %|| 'ENTHALPY'
% Enthalpy values can be in vaporization/condensation region    

% Initialize H and T values
    Hhi = T_or_H_hi;
    Hci = T_or_H_ci;
    Hhi_Ptr = libpointer('doublePtr',Hhi);
    Hci_Ptr = libpointer('doublePtr',Hci);
    Thi = zeros(var_length,1); %initialized to 0 temporarily
    Tci = zeros(var_length,1); %initialized to 0 temporarily
    Thi_Ptr = libpointer('doublePtr',Thi);
    Tci_Ptr = libpointer('doublePtr',Tci);

% Update values located at Thi and Tci pointer locations
    calllib('coolprop','AbstractState_update_and_1_out', fluid_h, hp_ind, Hhi_Ptr, p_h_Ptr, var_length, t_ind, Thi_Ptr, ierr,herr,buffer_size);
    calllib('coolprop','AbstractState_update_and_1_out', fluid_c, hp_ind, Hci_Ptr, p_c_Ptr, var_length, t_ind, Tci_Ptr, ierr,herr,buffer_size);
% Set Thi & Tci to values stored at pointer location
    Thi = get(Thi_Ptr,'Value');
    Tci = get(Tci_Ptr,'Value');
    
else
    error('Please specify whether input uses variables of Temperature or Enthalpy.')
end
    

% Determine vaporization temperature for fluids at respective pressures.
% Initialize pointer locations for relevant variables
Q_Ptr = libpointer('doublePtr', zeros(var_length,1)); % Q is specifically set to 0
Tvh = zeros(var_length,1); % initialized to 0 temporarily
Tvc = zeros(var_length,1); % initialized to 0 temporarily
Tvh_Ptr = libpointer('doublePtr', Tvh);
Tvc_Ptr = libpointer('doublePtr', Tvc);

% Update values located at the Tvc & Tvh pointer locations for vaporization
% conditions
calllib('coolprop','AbstractState_update_and_1_out', fluid_h, pq_ind, p_h_Ptr, Q_Ptr, var_length, t_ind, Tvh_Ptr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', fluid_c, pq_ind, p_c_Ptr, Q_Ptr, var_length, t_ind, Tvc_Ptr, ierr,herr,buffer_size);

% Set Tvh & Tvc to values stored at pointer location
Tvh = get(Tvh_Ptr,'Value');
Tvc = get(Tvc_Ptr,'Value');

% Determine vaporization enthalpies for fluids at respective pressures
% Min enthalpies of vaporization region for hot and cool fluids (Q = 0)
Hvh_min = zeros(var_length,1); % initialized to 0 temporarily
Hvc_min = zeros(var_length,1); % initialized to 0 temporarily
Hvh_min_Ptr = libpointer('doublePtr',Hvh_min);
Hvc_min_Ptr = libpointer('doublePtr',Hvc_min);
    
% Update values located at Hvh_max and Hvc_max pointer locations
calllib('coolprop','AbstractState_update_and_1_out', fluid_h, pq_ind, p_h_Ptr, Q_Ptr, var_length, h_ind, Hvh_min_Ptr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', fluid_c, pq_ind, p_c_Ptr, Q_Ptr, var_length, h_ind, Hvc_min_Ptr, ierr,herr,buffer_size);
% Set Hvh_max & Hvc_max to values stored at pointer location
Hvh_min = get(Hvh_min_Ptr,'Value');
Hvc_min = get(Hvc_min_Ptr,'Value');

% Max enthalpies of vaporization region for hot and cool fluid (Q = 1)
Q_Ptr = libpointer('doublePtr', ones(var_length,1)); % Q is specifically set to 1
Hvh_max = zeros(var_length,1); % initialized to 0 temporarily
Hvc_max = zeros(var_length,1); % initialized to 0 temporarily
Hvh_max_Ptr = libpointer('doublePtr',Hvh_max);
Hvc_max_Ptr = libpointer('doublePtr',Hvc_max);
    
% Update values located at Hvh_max and Hvc_max pointer locations
calllib('coolprop','AbstractState_update_and_1_out', fluid_h, pq_ind, p_h_Ptr, Q_Ptr, var_length, h_ind, Hvh_max_Ptr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', fluid_c, pq_ind, p_c_Ptr, Q_Ptr, var_length, h_ind, Hvc_max_Ptr, ierr,herr,buffer_size);
% Set Hvh_max & Hvc_max to values stored at pointer location
Hvh_max = get(Hvh_max_Ptr,'Value');
Hvc_max = get(Hvc_max_Ptr,'Value');


% Determine the specific heats of the incoming fluids using coolprop
% Initialize pointer locations for relevant variables
Ch = zeros(var_length,1); % initialized to 0 temporarily
Cc = zeros(var_length,1); % initialized to 0 temporarily
Ch_Ptr = libpointer('doublePtr', Ch); 
Cc_Ptr = libpointer('doublePtr', Cc); 

% Update values located at the Ch & Cc pointer locations for input
% conditions
calllib('coolprop','AbstractState_update_and_1_out', fluid_h, hp_ind, Hhi_Ptr, p_h_Ptr, var_length, c_ind, Ch_Ptr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', fluid_c, hp_ind, Hci_Ptr, p_c_Ptr, var_length, c_ind, Cc_Ptr, ierr,herr,buffer_size);

% Set Ch & Cc to values stored at pointer location
Ch = get(Ch_Ptr,'Value');   % J/(kg K)
Cc = get(Cc_Ptr,'Value');   % J/(kg K)

% If the fluid is withing the vaporization range then the coolprop library 
% will return an inaccurate values. If that happens then the specific heat 
% will be set to inf so its temperature will not change even though heat is
% exchanged.
Ch(Hhi <= Hvh_max & Hhi > Hvh_min) = inf;    % J/(kg K)
Cc(Hci < Hvc_max & Hci >= Hvc_min) = inf;    % J/(kg K)

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

% unitA is the size of heat exchanger area(m2) each iteration of the script
% uses for calculating NTU
% unitA = 1;
unitA = A .* ones(var_length,1);

% Calculate the number of transfer units of the heat exchanger
if any(A < unitA)
    unitA = A .* ones(var_length,1);
%     NTU = U * A / Cmin
end

NTU = U * unitA ./ Cmin;


% HeatExType = 'counter';
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
% if Cr == 0
%     effectiveness = 1 - exp(-NTU);
% end
% effectiveness;


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

% Determine Hho and Hco values based off of Tho and Tco values
% Initialize Hho and Hco values and pointers and Tho & Tco pointers
Tho_Ptr = libpointer('doublePtr',Tho);
Tco_Ptr = libpointer('doublePtr',Tco);
Hho = zeros(var_length,1); % initialized to 0 temporarily
Hco = zeros(var_length,1); % initialized to 0 temporarily
Hho_Ptr = libpointer('doublePtr',Hho);
Hco_Ptr = libpointer('doublePtr',Hco);
    
% Update values located at Hho and Hco pointer locations
calllib('coolprop','AbstractState_update_and_1_out', fluid_h, pt_ind, p_h_Ptr, Tho_Ptr, var_length, h_ind, Hho_Ptr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', fluid_c, pt_ind, p_c_Ptr, Tco_Ptr, var_length, h_ind, Hco_Ptr, ierr,herr,buffer_size);
% Set Hho & Hco to values stored at pointer location
Hho = get(Hho_Ptr,'Value');
Hco = get(Hco_Ptr,'Value');

% The coolprop will fail to update if Tho = vaporization temperature. In
% that case, those values will remain zero because it was initialized to 
% zero. For those values Ho = Hi - P/m_dot
Hho(~Hho) =  Hhi(~Hho) - P(~Hho)./m_dot_h(~Hho);

% The coolprop will fail to update if Tco = vaporization temperature. In
% that case, those values will remain zero because it was initialized to 
% zero. For those values Ho = Hi + P/m_dot
Hco(~Hco) =  Hci(~Hco) + P(~Hco)./m_dot_c(~Hco);


% check if hot fluid begins to condense 
if any(Hho <= Hvh_max & Hhi > Hvh_max)
% if Tho <= Tvh && Thi >= Tvh
%     disp('~~~~~~~~~~~~~~~~~~~Condensation start~~~~~~~~~~~~~~~~~~~~~')

%     Get the indecies where condensation starts
    cond_start_ind = (Hho <= Hvh_max & Hhi > Hvh_max);
    
%     Set the H output of the appropriate values to maximum value in the
%     vaporization region (Q = 1)
    Hho(cond_start_ind) = Hvh_max(cond_start_ind);
    Tho(cond_start_ind) = Tvh(cond_start_ind);
    
%     Calculate power needed to bring hot fluid initial state to
%     precondense for only relevant values
    P_init_to_precondense = (Hhi(cond_start_ind) - Hho(cond_start_ind)).*m_dot_h(cond_start_ind);

%     Calculate change in temperature for relavent cool fluid values
    deltaTc = P_init_to_precondense./Cc_rate(cond_start_ind);
%     Update relavent cool fluid output temperatures and get new pointer
    Tco(cond_start_ind) = Tci(cond_start_ind) + deltaTc;
    Tco_Ptr = libpointer('doublePtr',Tco);
    
%     Update H values for new cool output temperatures
    calllib('coolprop','AbstractState_update_and_1_out', fluid_c, pt_ind, p_c_Ptr, Tco_Ptr, var_length, h_ind, Hco_Ptr, ierr,herr,buffer_size);
    Hco = get(Hco_Ptr,'Value');

%     Calculate ratio of power to precondense state relative to initally
%     calculated power in order to correct the area of relavent values
    PowerFrac_Condense_over_NoCondense = P_init_to_precondense./P(cond_start_ind);
    unitA(cond_start_ind) = unitA(cond_start_ind)*PowerFrac_Condense_over_NoCondense;
%     Update relevant power values
    P(cond_start_ind) = P_init_to_precondense;
    
% check if hot fluid finishes condensing    
elseif any(Hho <= Hvh_min & Hhi > Hvh_min)
% elseif Tho <= Tvh && Thi >= Tvh
%     disp('~~~~~~~~~~~~~~~~~~~~~~Condensation end~~~~~~~~~~~~~~~~~~~~~~~~~')

%     Get the indecies where condensation starts
    cond_end_ind = (Hho <= Hvh_min & Hhi > Hvh_min);
    
%     Set the H output of the appropriate values to minimum value in the
%     vaporization region (Q = 0)
    Hho(cond_end_ind) = Hvh_min(cond_end_ind);
    Tho(cond_end_ind) = Tvh(cond_end_ind) - 1e-4;
    
%     Calculate power needed to bring hot fluid initial state to
%     postcondense for only relevant values
    P_init_to_postcondense = (Hhi(cond_end_ind) - Hho(cond_end_ind)).*m_dot_h(cond_end_ind);

%     Calculate change in temperature for relavent cool fluid values
    deltaTc = P_init_to_postcondense./Cc_rate(cond_end_ind);
%     Update relavent cool fluid output temperatures and get new pointer
    Tco(cond_end_ind) = Tci(cond_end_ind) + deltaTc;
    Tco_Ptr = libpointer('doublePtr',Tco);

%     Update H values for new cool output temperatures
    calllib('coolprop','AbstractState_update_and_1_out', fluid_c, pt_ind, p_c_Ptr, Tco_Ptr, var_length, h_ind, Hco_Ptr, ierr,herr,buffer_size);
    Hco = get(Hco_Ptr,'Value');

%     Calculate ratio of power to postcondense state relative to initally
%     calculated power in order to correct the area of relavent values
    PowerFrac_Condense_over_NoCondense = P_init_to_postcondense./P(cond_end_ind);
    unitA(cond_end_ind) = unitA(cond_end_ind)*PowerFrac_Condense_over_NoCondense;
%     Update relevant power values
    P(cond_end_ind) = P_init_to_postcondense;
 

% else check if cool fluid begins to vaporize
elseif any(Hco >= Hvc_min & Hci < Hvc_min)
% elseif Tco >= Tvc && Tci <= Tvc
%    disp('~~~~~~~~~~~~~~~~~~~Vaporization start~~~~~~~~~~~~~~~~~~~~~~~~')

%     Get the indecies where vaporization starts
    vap_start_ind = (Hco >= Hvc_min & Hci < Hvc_min);

%     Set the H output of the appropriate values to minimum value in the
%     vaporization region (Q = 0)
    Hco(vap_start_ind) = Hvc_min(vap_start_ind);
    Tco(vap_start_ind) = Tvc(vap_start_ind);

%     Calculate power needed to bring cool fluid initial state to
%     prevaporization for only relevant values
    P_init_to_prevap = (Hco(vap_start_ind) - Hci(vap_start_ind)).*m_dot_c(vap_start_ind);

%     Calculate change in temperature for relavent hot fluid values
    deltaTh = P_init_to_prevap./Ch_rate(vap_start_ind);
%     Update relavent hot fluid output temperatures and get new pointer
    Tho(vap_start_ind) = Thi(vap_start_ind) - deltaTh;
    Tho_Ptr = libpointer('doublePtr',Tho);

%     Update H values for new hot output temperatures
    calllib('coolprop','AbstractState_update_and_1_out', fluid_h, pt_ind, p_h_Ptr, Tho_Ptr, var_length, h_ind, Hho_Ptr, ierr,herr,buffer_size);
    Hho = get(Hho_Ptr,'Value');

%     Calculate ratio of power to prevaporize state relative to initally
%     calculated power in order to correct the area of relavent values
    PowerFrac_Vaporize_over_NoVaporize = P_init_to_prevap./P(vap_start_ind);
    unitA(vap_start_ind) = unitA(vap_start_ind).*PowerFrac_Vaporize_over_NoVaporize;
%     Update relevant power values
    P(vap_start_ind) = P_init_to_prevap;

    
% check if cool fluid finishes vaporizing
elseif any(Hco >= Hvc_max & Hci < Hvc_max)
% elseif Tho <= Tvh && Thi >= Tvh
%     disp('~~~~~~~~~~~~~~~~Vaporization end~~~~~~~~~~~~~~~~~~~~')

%     Get the indecies where vaporization ends
    vap_end_ind = (Hco >= Hvc_max & Hci < Hvc_max);

%     Set the H output of the appropriate values to maximum value in the
%     vaporization region (Q = 1)
    Hco(vap_end_ind) = Hvc_max(vap_end_ind);
    Tco(vap_end_ind) = Tvc(vap_end_ind) + 1e-4;

%     Calculate power needed to bring cool fluid initial state to
%     postvaporization for only relevant values
    P_init_to_postvap = (Hco(vap_end_ind) - Hci(vap_end_ind)).*m_dot_c(vap_end_ind);

%     Calculate change in temperature for relavent hot fluid values
    deltaTh = P_init_to_postvap./Ch_rate(vap_end_ind);
%     Update relavent hot fluid output temperatures and get new pointer
    Tho(vap_end_ind) = Thi(vap_end_ind) - deltaTh;
    Tho_Ptr = libpointer('doublePtr',Tho);

%     Update H values for new hot output temperatures
    calllib('coolprop','AbstractState_update_and_1_out', fluid_h, pt_ind, p_h_Ptr, Tho_Ptr, var_length, h_ind, Hho_Ptr, ierr,herr,buffer_size);
    Hho = get(Hho_Ptr,'Value');

%     Calculate ratio of power to prevaporize state relative to initally
%     calculated power in order to correct the area of relavent values
    PowerFrac_Vaporize_over_NoVaporize = P_init_to_postvap./P(vap_end_ind);
    unitA(vap_end_ind) = unitA(vap_end_ind).*PowerFrac_Vaporize_over_NoVaporize;
%     Update relevant power values
    P(vap_end_ind) = P_init_to_postvap;

else
%     disp('~~~~~~~~~~~~~~~~Maintain State~~~~~~~~~~~~~~~~~~~~')
end
%     Hhi
%     Hho
%     Thi
%     Tho
%     Hci
%     Hco
%     Tci
%     Tco
%     P

% double check if hot fluid would condense and cool fluid would vaporize
if any((Tho < Tvh & Thi > Tvh) & (Tco > Tvc & Tci < Tvc))
    error('This function cannot accurately handle both fluids changing phase.')
end

if any(A > unitA)

%     P_next = zeros(size(T_or_H_hi));
%     disp(['Power Before: ' num2str(P) ' W'])
    [ P_next, Hho, Hco ] = HeatExNTU_3_low( Hho, m_dot_h, p_h, ...
        Hco, m_dot_c, p_c, fluid_h, fluid_c, HeatExType, U, A - unitA, 'H' );
%     A
%     unitA
% %     A - unitA;
%     P_next
    P = P + P_next;
%     disp(['Power After: ' num2str(P) ' W'])

    Hho_Ptr = libpointer('doublePtr',Hho);
    Hco_Ptr = libpointer('doublePtr',Hco);
% Update values located at Tho and Tco pointer locations
    calllib('coolprop','AbstractState_update_and_1_out', fluid_h, hp_ind, Hho_Ptr, p_h_Ptr, var_length, t_ind, Tho_Ptr, ierr,herr,buffer_size);
    calllib('coolprop','AbstractState_update_and_1_out', fluid_c, hp_ind, Hco_Ptr, p_c_Ptr, var_length, t_ind, Tco_Ptr, ierr,herr,buffer_size);
% Set Thi & Tci to values stored at pointer location
    Tho = get(Tho_Ptr,'Value');
    Tco = get(Tco_Ptr,'Value');

% else
%     P
end


% Free the hot and cool fluid handles
% calllib('coolprop','AbstractState_free', fluid_h, ierr,herr,buffer_size);
% calllib('coolprop','AbstractState_free', fluid_c, ierr,herr,buffer_size);

% if upper(T_or_H) == 'T' %|| 'TEMPERATURE'
if T_or_H == 'T' %|| 'TEMPERATURE'
    T_or_H_ho = Tho;
    T_or_H_co = Tco;
% elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
elseif T_or_H == 'H' %|| 'ENTHALPY'
    T_or_H_ho = Hho;
    T_or_H_co = Hco;
end
end

