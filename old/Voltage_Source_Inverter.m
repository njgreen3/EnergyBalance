function [ Ploss, Qloss ] = Voltage_Source_Inverter(Pload, Qload, Vdc, VacOut, VacIn, switchFreq, ...
                                VthreshTran, VthreshDiode, RcondTran, RcondDiode,...
                                transistorNum, diodeNum, duty )
                            
%INVERTER Summary of this function goes here
%   This function simulates a three phase voltage source inverter for an
%   energy balance microgrid model.
%   The losses in a switching based inverter are generally broken up into
%   losses during conduction through the switching device and losses during
%   the switch transition (turning on and off). Switching devices include
%   transistors and diodes.

% The calculations in this function are based primarily on those found in
% the article Comparison of Power Losses, Current and Voltage Stresses of
% Semiconductors in Voltage Source Transformerless Multilevel Inverters
% by Perantzakis et. al. for the journal IET Power Electronics

end

