function [ P, T_or_H_ho, T_or_H_co ] = HeatExNTU_3( T_or_H_hi, m_dot_h, p_h, fluid_h, ...
    T_or_H_ci, m_dot_c, p_c, fluid_c, HeatExType, U, A, T_or_H )
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


% Ensure temperature/enthalpy input variables are delt with appropriately
if upper(T_or_H) == 'T' %|| 'TEMPERATURE'
% Assume given temperature values are not at vaporization/condensation
    Thi = T_or_H_hi;
    Tci = T_or_H_ci;
    Hhi = CoolProp.PropsSI('H','P',p_h,'T',Thi,fluid_h);
    Hci = CoolProp.PropsSI('H','P',p_c,'T',Tci,fluid_c);
elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
% Enthalpy values can be in vaporization/condensation region    
    Hhi = T_or_H_hi;
    Hci = T_or_H_ci;
    Thi = CoolProp.PropsSI('T','P',p_h,'H',Hhi,fluid_h);
    Tci = CoolProp.PropsSI('T','P',p_c,'H',Hci,fluid_c);
else
    error('Please specify whether input uses variables of Temperature or Enthalpy.')
end
    

% Determine vaporization temperature for fluids at respective pressures
Tvh = CoolProp.PropsSI('T','P',p_h,'Q',0,fluid_h);
Tvc = CoolProp.PropsSI('T','P',p_c,'Q',0,fluid_c);

% Determine vaporization enthalpies for fluids at respective pressures
% Max and min enthalpies of vaporization region for hot fluid
Hvh_max = CoolProp.PropsSI('H','P',p_h,'Q',1,fluid_h);
Hvh_min = CoolProp.PropsSI('H','P',p_h,'Q',0,fluid_h);
% Max and min enthalpies of vaporization region for cool fluid
Hvc_max = CoolProp.PropsSI('H','P',p_c,'Q',1,fluid_c);
Hvc_min = CoolProp.PropsSI('H','P',p_c,'Q',0,fluid_c);

% Determine the specific heats of the incoming fluids using coolprop
% If the fluid is withing 1e-4 units of the vaporization point then the 
% command will fail. If that happens then the specific heat will be set to 
% inf so its temperature will not change even though heat is exchanged.
try
    Ch = CoolProp.PropsSI('C', 'H', Hhi, 'P', p_h, fluid_h);    % J/(kg K)
catch %ME
    Ch = inf;
end
try
    Cc = CoolProp.PropsSI('C', 'H', Hci, 'P', p_c, fluid_c);    % J/(kg K)
catch %ME
%     w = warning('query','last')
    Cc = inf;
end
Ch(Hhi <= Hvh_max & Hhi > Hvh_min) = inf;    % J/(kg K)
Cc(Hci < Hvc_max & Hci >= Hvc_min) = inf;    % J/(kg K)

% Calculate the heat capacity rates of the hot and cool fluids
Ch_rate = Ch * m_dot_h;     % W/K
Cc_rate = Cc * m_dot_c;     % W/K

% Determine the maximum and minimum of the two heat capacity rates
Cmax = max(Ch_rate, Cc_rate);   % W/K
Cmin = min(Ch_rate, Cc_rate);   % W/K

% Calcluate the heat capacity rate ratio 
Cr = Cmin/Cmax;

% unitA is the size of heat exchanger area(m2) each iteration of the script
% uses for calculating NTU
% unitA = 1;
unitA = A;

% Calculate the number of transfer units of the heat exchanger
if A < unitA
    unitA = A;
%     NTU = U * A / Cmin
end

NTU = U * unitA / Cmin;


% Calculate effectiveness based on type of heat exchanger
switch HeatExType
    case 'counter'  % If the heat exchanger is counter-flow
        effectiveness = (1 - exp(-NTU.*(1 - Cr))) ./ ...
                        (1 - Cr.*exp(-NTU.*(1 - Cr)));
    case 'parallel' % If the heat exchanger is parallel-flow
        effectiveness = (1 - exp(-NTU.*(1 + Cr))) ./ ...
                        (1 + Cr);
    otherwise       % If not exchanger specified
        warning('No heat exchanger type specified, therefore no heat is exchanged.')
        effectiveness = 0;
end        

% However if the Cr is 0 (Cmax is inf) then effectiveness is calcluated
% regardless of the heat exchanger flow.
if Cr == 0
    effectiveness = 1 - exp(-NTU);
end
effectiveness;
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

% The cool prop function may fail if Tho = vaporization temperature
try
%     disp('try')
    Hho = CoolProp.PropsSI('H','P',p_h,'T',Tho,fluid_h);
catch %ME
%     disp('fail')
    Hho = Hhi - P./m_dot_h;
end
% The cool prop function may fail if Tco = vaporization temperature
try
    Hco = CoolProp.PropsSI('H','P',p_c,'T',Tco,fluid_c);
catch %ME
    Hco = Hci + P./m_dot_c;
end
    
% check if hot fluid begins to condense 
if Hho <= Hvh_max && Hhi > Hvh_max
% if Tho <= Tvh && Thi >= Tvh
%     disp('~~~~~~~~~~~~~~~~~~~Condensation start~~~~~~~~~~~~~~~~~~~~~')
    Hho = Hvh_max;
    Tho = Tvh - 1e-4;
        
    P_init_to_precondense = (Hhi - Hho).*m_dot_h;
       
    deltaTc = P_init_to_precondense./Cc_rate;
    Tco = Tci + deltaTc;
    Hco = CoolProp.PropsSI('H','P',p_c,'T',Tco,fluid_c);
    
    PowerFrac_Condense_over_NoCondense = P_init_to_precondense/P;
    
    unitA = unitA*PowerFrac_Condense_over_NoCondense;
   
    P = P_init_to_precondense;
    
% check if hot fluid finishes condensing    
elseif Hho <= Hvh_min && Hhi > Hvh_min
% elseif Tho <= Tvh && Thi >= Tvh
%     disp('~~~~~~~~~~~~~~~~~~~~~~Condensation end~~~~~~~~~~~~~~~~~~~~~~~~~')
    Hho = Hvh_min;
    Tho = Tvh;
    
    P_init_to_postcondense = (Hhi - Hho).*m_dot_h;
    
    deltaTc = P_init_to_postcondense./Cc_rate;
    Tco = Tci + deltaTc;
    Hco = CoolProp.PropsSI('H','P',p_c,'T',Tco,fluid_c);
    
    PowerFrac_Condense_over_NoCondense = P_init_to_postcondense/P;
    
    unitA = unitA*PowerFrac_Condense_over_NoCondense;
   
    P = P_init_to_postcondense;
    

% else check if cool fluid begins to vaporize
elseif Hco >= Hvc_min && Hci < Hvc_min
% elseif Tco >= Tvc && Tci <= Tvc
%    disp('~~~~~~~~~~~~~~~~~~~Vaporization start~~~~~~~~~~~~~~~~~~~~~~~~')

    Hco = Hvc_min;
    Tco = Tvc;
        
    P_init_to_prevap = (Hco - Hci).*m_dot_c;
       
    deltaTh = P_init_to_prevap./Ch_rate;
    Tho = Thi - deltaTh;
    Hho = CoolProp.PropsSI('H','P',p_c,'T',Tho,fluid_h);
    
    PowerFrac_Vaporize_over_NoVaporize = P_init_to_prevap/P;
    
    unitA = unitA*PowerFrac_Vaporize_over_NoVaporize;
   
    P = P_init_to_prevap;
    
% check if cool fluid finishes vaporizing
elseif Hco >= Hvc_max && Hci < Hvc_max
% elseif Tho <= Tvh && Thi >= Tvh
%     disp('~~~~~~~~~~~~~~~~Vaporization end~~~~~~~~~~~~~~~~~~~~')
    Hco = Hvc_max;
    Tco = Tvc + 1e-4;
    
    P_init_to_postvap = (Hco - Hci).*m_dot_c;
    
    deltaTh = P_init_to_postvap./Ch_rate;
    Tho = Thi - deltaTh;
    Hho = CoolProp.PropsSI('H','P',p_h,'T',Tho,fluid_h);
    
    PowerFrac_Vaporize_over_NoVaporize = P_init_to_postvap/P;
    
    unitA = unitA*PowerFrac_Vaporize_over_NoVaporize;
   
    P = P_init_to_postvap;
    
end

% double check if hot fluid would condense and cool fluid would vaporize
if (Tho < Tvh && Thi > Tvh) && (Tco > Tvc && Tci < Tvc)
    error('This function cannot accurately handle both fluids changing phase.')
end

if A > unitA

%     disp(['Power Before: ' num2str(P) ' W'])
    [ P_next, Hho, Hco ] = HeatExNTU_3( Hho, m_dot_h, p_h, fluid_h, ...
        Hco, m_dot_c, p_c, fluid_c, HeatExType, U, A - unitA, 'H' );
    A;
    unitA;
%     A - unitA;
    P;
    P = P + P_next;
%     disp(['Power After: ' num2str(P) ' W'])
    
    Tho = CoolProp.PropsSI('T','P',p_c,'H',Hho,fluid_h);
    Tco = CoolProp.PropsSI('T','P',p_c,'H',Hco,fluid_c);
    
    
% else
%     P
end



if upper(T_or_H) == 'T' %|| 'TEMPERATURE'
    T_or_H_ho = Tho;
    T_or_H_co = Tco;
elseif upper(T_or_H) == 'H' %|| 'ENTHALPY'
    T_or_H_ho = Hho;
    T_or_H_co = Hco;
end
end

