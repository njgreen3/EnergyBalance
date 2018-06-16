function [ P, Tho, Tco ] = HeatExNTU_2( Thi, m_dot_h, p_h, fluid_h, ...
    Tci, m_dot_c, p_c, fluid_c, HeatExType, U, A )
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

% Determine vaporization temperature for fluids at respective pressures
Tvh = CoolProp.PropsSI('T','P',p_h,'Q',0,fluid_h);
Tvc = CoolProp.PropsSI('T','P',p_c,'Q',0,fluid_c);

% Determine the specific heats of the incoming fluids using coolprop
% If the fluid is withing 1e-4 units of the vaporization point then the 
% command will fail. If that happens then the specific heat will be set to 
% inf so its temperature will not change even though heat is exchanged.
try
    Ch = CoolProp.PropsSI('C', 'T', Thi, 'P', p_h, fluid_h);    % J/(kg K)
catch ME
    Ch = inf;
end
try
    Cc = CoolProp.PropsSI('C', 'T', Tci, 'P', p_c, fluid_c);    % J/(kg K)
catch ME
    Cc = inf;
end

% Calculate the heat capacity rates of the hot and cool fluids
Ch_rate = Ch * m_dot_h;     % W/K
Cc_rate = Cc * m_dot_c;     % W/K

% Determine the maximum and minimum of the two heat capacity rates
Cmax = max(Ch_rate, Cc_rate);   % W/K
Cmin = min(Ch_rate, Cc_rate);   % W/K

% Calcluate the heat capacity rate ratio 
Cr = Cmin/Cmax;

% Calculate the number of transfer units of the heat exchanger
NTU = U * A / Cmin;

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
if Tho < Tvh && Thi > Tvh
%    disp('Condensation')
%  Determine the enthalpy extracted from the hot fluid to condense it from 
%  starting temeprature. A slight 1e-4 K temperature change beyond
%  condensation is subtracted because coolprop cannot handle fluids within 
%  1e-4 K of the vaporization/condensation curve.
   Hh_init = CoolProp.PropsSI('H', 'T', Thi, 'P', p_h, fluid_h);
   Hh_cond = CoolProp.PropsSI('H', 'T', Tvh - 1e-4, 'P', p_h, fluid_h);
      
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
   [ P_deltaT, Tho, Tco ] = HeatExNTU_2( Tvh - 1e-4, m_dot_h, p_h, fluid_h, ...
       Tc_cond, m_dot_c, p_c, fluid_c, HeatExType, U, A/2 );
   P = P_cond + P_deltaT;

% else check if cool fluid would vaporize
elseif Tco > Tvc && Tci < Tvc
%    disp('Vaporization')
%  Determine the enthalpy needed to vaporize the cool fluid from 
%  starting temeprature. A slight 1e-4 K temperature change beyond
%  vaporization is added because coolprop cannot handle fluids within 
%  1e-4 K of the vaporization/condensation curve.
   Hc_init = CoolProp.PropsSI('H', 'T', Tci, 'P', p_c, fluid_c);
   Hc_vap = CoolProp.PropsSI('H', 'T', Tvc + 1e-4, 'P', p_c, fluid_c);
      
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
   [ P_deltaT, Tho, Tco ] = HeatExNTU_2( Th_vap, m_dot_h, p_h, fluid_h, ...
       Tvc + 1e-4, m_dot_c, p_c, fluid_c, HeatExType, U, A/2 );
   P = P_vap + P_deltaT;
end

% check if hot fluid would condense and cool fluid would vaporize
if (Tho < Tvh && Thi > Tvh) && (Tco > Tvc && Tci < Tvc)
    error('This function cannot accurately handle both fluids changing phase.')
end

end

