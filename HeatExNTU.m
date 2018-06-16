function [ P, Tho, Tco ] = HeatExNTU( Thi, m_dot_h, Ch, Tci, m_dot_c, Cc, ...
    HeatExType, U, A )
%HeatExNTU Summary 
%   This function calculates heat flow and outlet temperatures of a heat exhanger
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
% Note: The effectiveness is separte from efficiency. The calculations
%   assumes an efficiency of 1 where all heat is transferred from the hot
%   fluid to the cool fluid or remains within the the fluid. This also
%   assumes no fluid is undergoing a change in state, which may not be
%   sufficient for my final model.
% 
%   This funtion is unable to account for or detect phase changes. It 
%   should not be used for calculating heat transfer rates or temperature
%   changes in cases of evaporation or condensation. 
%
% Inputs:
%   Thi         Hot fluid inlet temperature in Celsius or Kelvin
%   m_dot_h     Hot fluid mass flow rate in kg/s
%   Ch          Hot fluid specfic heat in J/(kg K)
%   Tci         Cool fluid inlet temperature in Celsius or Kelvin
%   m_dot_c     Cool fluid mass flow rate in kg/s
%   Cc          Cool fluid specfic heat in J/(kg K)
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
%   Tho         Hot fluid outlet temperature in input temp units (C or K)
%   Tco         Cool fluid outlet temperature in input tempt units (C or K)

% Calculate the heat capacity rates of the hot and cool fluids
Ch_rate = Ch * m_dot_h;
Cc_rate = Cc * m_dot_c;

% Determine the maximum and minimum of the two heat capacity rates
Cmax = max(Ch_rate, Cc_rate);
Cmin = min(Ch_rate, Cc_rate);

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
        warning('No heat exchanger type specified, therefore no heat is exchanged.')
        effectiveness = 0;
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

end

