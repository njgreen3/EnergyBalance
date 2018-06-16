addpath('C:\Users\njgreen3\Documents\MATLAB\CoolPropStuff')
clear
warning('off','SWIG:RuntimeError')

% %%
% This script is a partial test of the ORC_Energy_Balance function.
% The script cycles through many different pressure values (both high and
% low) in order to determine pressure pairs which result in pump output
% temperature which is approximately equal to the initial temperature. The
% script records the pressure values into the arrays low_p and high_p.


% Define fluids
workingFluid = 'R245fa';
sourceFluid = 'Water';
sinkFluid = 'Water';

% Mass flow rates in kg/s % 12 kg/s is just shy of 200 gal/min of water
% m_dot_working = 1;
m_dot_source = 6;
m_dot_sink = 6;

% Temperatures of source, sink, and working fluid in K
T_source = 95 + 273;
T_sink = 5 + 273;
% T_working_init = 20 + 273;

% one atmosphere of pressure in pa for the pressure of the source and sink
atm = 101325;

% High and low pressure of working fluid in Pa
% pressure_h = 6.95e5; %7.0e5;
% pressure_l = 1.25e5; %1.1e5;
pressure_source = atm;
pressure_sink = atm;

% Set pump and expander isentropic efficiencies
eff_pump = .7;
eff_expander = .78;

% heat exchanger parameters for evaporator
heatExTypeVap = 'counter';
heatExAVap = 9.3; %268.5
heatExUVap = 1500;
    
% heat exchanger parameters for condensor
heatExTypeCond = 'counter';
heatExACond = 18.6; %268.5
heatExUCond = 1400;

% initialize count number of pressure combinations where T_working_init 
% approximately equals T_pump_out
n = 0;

% the error value between T_working_init and T_pump_out such that they be
% considered approximately equal
err_T = 0.5;

% maximum pressure differential across expander
max_pres_diff = 10e5;

% cycle through high pressures
for pressure_h = 2.0e5:10e4:10.0e5
%   cycle through low pressures ensuring low pressure is always lower than high pressure
%     disp(['high pressure ' num2str(pressure_h) ' pa'])
    for pressure_l = 1.5e5:5e4:pressure_h-5e4
%         disp(['low pressure ' num2str(pressure_l) ' pa'])
%       Try ORC function. Likely errors include trying to use Coolprop to
%       analyze a fluid directly on a vaporization curve because it 
%       cannot determine phase with given info

        if pressure_h - pressure_l > max_pres_diff
            continue
        end

        for m_dot_working = (1:10)%*0.5
%             disp(['flow rate ' num2str(m_dot_working) ' kg/s'])
        for T_working_init = (7:2:27) + 273
%             disp(['temperature ' num2str(T_working_init) ' K'])
        try     
            [ P, T ] = ORC_Energy_Balance( ...
                workingFluid, sourceFluid, sinkFluid, m_dot_working, m_dot_source, m_dot_sink, ...
                T_source, T_sink, T_working_init, pressure_h, pressure_l, pressure_source, pressure_sink,...
                eff_pump, eff_expander, ...
                heatExTypeVap, heatExAVap, heatExUVap ,heatExTypeCond, heatExACond, heatExUCond);
        catch ME
%           If it does encounter an error, just continue to the next prssure value
%             ME
%             ME.stack.name
%             ME.stack.line

%             disp(['The function failed at ' num2str(pressure_l) ' and ' num2str(pressure_h) ' pa']);
            continue
%             break
        end
%         disp(n)
%         disp(num2str(T.pump_out))
%       Test if T_pump_out is within err of T_working_init
        if abs(T.working_init - T.pump_out) < err_T
            if abs(P.outMech) > abs(P.inMech)
%              Increase count
                n = n + 1;
                disp(num2str(n))
%               Record successful pressures into arrays
                low_p(n) = pressure_l;
                high_p(n) = pressure_h;
                P_out(n) = P.outMech;
                P_in(n) = -P.inMech;
                T_pump(n) = T.working_init;
                flow(n) = m_dot_working;
                heat_in(n) = -P.inHeat;
                heat_out(n) = P.outHeat;
                source_out(n) = T.source_out;
                sink_out(n) = T.sink_out;
            end
        end
        end
        end
    end
end
P_e_out = P_out*0.9;
P_e_in = P_in/0.9;
P_net = P_out - P_in;
P_e_net = P_e_out - P_e_in;

figure()
plot3(high_p-low_p,T_pump,P_e_out, '*')
grid
xlabel('\Delta Pressure (Pa)')
ylabel('Initial Temperature (K)')
zlabel('Power_e out (W)')

figure()
plot3(flow,T_pump,P_e_out, '*')
grid
xlabel('Working Fluid Flow (kg/s)')
ylabel('Initial Temperature (K)')
zlabel('Power_e out (W)')


figure()
subplot(2,2,4)
plot3(high_p,low_p,P_e_out, '*')
title({'d. Electrical Power Out vs. High and Low';'Pressure Set Points'})
grid
xlabel('High Pressure (Pa)')
ylabel('Low Pressure (Pa)')
zlabel('Power_e out (W)')


% figure()
subplot(2,2,2)
plot3(high_p-low_p,flow,P_e_out, '*')
title({'b. Electrical Power Out vs. Working Fluid'; 'Pressure Drop and Flow Rate'})
grid
xlabel('\Delta Pressure (Pa)')
ylabel('Working Fluid Flow (kg/s)')
zlabel('Power_e out (W)')

% figure()
subplot(2,2,3)
plot3(P_e_out,P_e_in,P_e_net, '*')
title({'c. Net Electrical Power vs. Electrical Power Out'; 'and Electrical Power In'})
grid
xlabel('Power_e out (W)')
ylabel('Working Fluid Pump Power_e (W)')
zlabel('Net Power_e (W)')

% figure()
subplot(2,2,1)
plot3(heat_out,heat_in,P_out, '*')
title({'a. Mechanical Power Out vs. Heat Flow'; ' Into and Out of the System'})
grid
xlabel('Heat expelled (W)')
ylabel('Heat in (W)')
zlabel('Power_m out (W)')

% figure()
% plot3(P_in,T_pump,P_net, '*')
% grid
% xlabel('Pump Power in (W)')
% ylabel('Initial Temperature (K)')
% zlabel('Net Power (W)')

warning('on','SWIG:RuntimeError')