function [efficiency] = InverterEfficiency(operatingPower, ratedPower)
%INVERTEREFFICIENCY 
%   This function takes the operating power and rated power of an inverter
%   as inputs and returns the operating efficiency. The efficiency is
%   estimated based off of the AEG Power Solutions model PV.500-UL. The
%   data are obtained from the gosolarcalifornia.ca.gov database of solar
%   inverters, specifically http://www.gosolarcalifornia.ca.gov/equipment/inverter_tests/summaries/AEG%20PV%20500-UL%20with%20Transformer.pdf
%   The efficiency values also vary with the DC voltage setpoint of the 
%   input with, datapoints at minimum, nominal, and maximum voltages. This 
%   script will assume nomimal voltage.

if nargin == 0
    operatingPower = 4e4;
    ratedPower = 5e4;
end

testPower = [ 0.1 0.2 0.3 0.5 0.75 1 ];
testEff = [ 0.937 0.963 0.971 0.975 0.976 0.974 ];

relPower = operatingPower/ratedPower;

efficiency = interp1(testPower, testEff, relPower, 'pchip');

end

