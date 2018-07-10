function [ Pcond, Psw ] = VSI_semiconductor_loss(Pload, Qload, Vac, switchFreq, ...
                                Vthresh, Vblock, Rcond, turnOnTime, turnOffTime, ...
                                duty, modAmp, swTheta1, swTheta2, condTheta1, condTheta2 )
                            
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


if nargin == 0
   Pload = 5000;    % W
   Qload = 2000;    % VAR
   Vac = 480;   % V
   switchFreq = 50e3;   % Hz
   Vthresh = 4; % V
   Vblock = 5;   % V
   Rcond = 1e-3; % ohm
   turnOnTime = 3e-6;   % s
   turnOffTime = 3e-6;  % s
   duty = .5;   
   modAmp = .8;
   swTheta1 = 0; %?
   swTheta2 = .1*pi; %?  
   condTheta1 = .1*pi;   %?
   condTheta2 = 1*pi;   %?

end

Imax = sqrt(2)*abs((Pload + 1i*Qload)/Vac);

pf = Pload/abs(Pload + 1i*Qload);
pf_angle = acos(pf);

iload = @(theta) Imax*modAmp*sin(theta - pf_angle);

Psw = 1/(2*pi) * Vblock*(turnOnTime + turnOffTime) * switchFreq/2 * ...
    integral(@(theta) abs(Imax*modAmp*sin(theta - pf_angle)), swTheta1, swTheta2);

Pcond = 1/(2*pi) * integral(@(theta) abs(Imax*modAmp*sin(theta - pf_angle)).*(Vthresh + Rcond*abs(Imax*modAmp*sin(theta - pf_angle)))*duty,...
    condTheta1,condTheta2);

end

