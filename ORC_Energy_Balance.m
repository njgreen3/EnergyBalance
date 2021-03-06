function [ ORC_Power, ORC_Temperature ] = ...
    ORC_Energy_Balance( workingFluid, sourceFluid, sinkFluid,...
    m_dot_working, m_dot_source, m_dot_sink, ...
    T_source, T_sink, T_working_init, pressure_h, pressure_l, ...
    pressure_source, pressure_sink, eff_pump, eff_expander, ...
    heatExTypeVap, heatExAVap, heatExUVap, heatExTypeCond, heatExACond, heatExUCond)
%ORC_ENERGY_BALANCE - Simulates an Organic rankine cycle.
%   The function makes use of the functions HeatExNTU_2 and IsentropicPower
%   in order calculate power values at different stages of the ORC cycle.
%   It is assumed no ambient heat is lost at any stage during the cycle.
%   It is assumed the evaporator and condensor are isobaric (no pressure
%       drop.
%   It is assumed the pump and expander are non-ideal isentropic processes.
%       This means power calculations are made assuming entropy remains 
%       constant during each process then the efficiency is applied and 
%       output power and temperature values are recalculated. 
% 
%   This function makes use of the CoolProp wrapper for fluid properties. 
%   The high level Matlab wrapper cannot handle vector inputs. There are 
%   ways to make low level calls interfacing with the DLL so many values 
%   can be determined more quickly, but that is not implented here.
%   More information about the Matlab wrapper can be found at
%     http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html
%   A list of parameters can be found at 
%     http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table
%   A list of fluids can be found at 
%     http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids

% Inputs:
%   workingFluid    Name of working fluid (string)
%   sourceFluid     Name of fluid used as heat source (string)
%   sinkFluid       Name of fluid used as heat sink (string)
%   m_dot_working   Working fluid mass flow rate in kg/s
%   m_dot_source    Source fluid mass flow rate in kg/s
%   m_dot_sink      Sink fluid mass flow rate in kg/s
%   T_source        Source fluid inlet temperature in K
%   T_sink          Sink fluid inlet temperature in K
%   T_working_init  Working fluid initial temperature in K
%   pressure_h      Working fluid pressure for high pressure region in pa
%   pressure_l      Working fluid pressure for low pressure region in pa
%   pressure_source Source fluid inlet and outlet pressure
%   pressure_sink   Sink fluid inlet and outlet pressure
%   eff_pump        Isentropic efficiency of the pump. Expected range: 0-1
%   eff_expander    Isentropic efficiency of the expander. Expected range: 0-1
%   HeatExTypeVap   String describing the Evaporator geometry. 
%                   Currently this variable only accepts 'counter' and 
%                   'parallel', but can be expanded to other flow types
%                   'counter' for counter flow
%                   'parallel' for parallel
%   HeatExUVap      Evaporator Overall heat transfer coefficient in W/(K m2)
%   HeatExAVap      Evaporator Heat transfer area in m2
%   HeatExTypeCond  String describing the Condensor geometry. 
%                   Currently this variable only accepts 'counter' and 
%                   'parallel', but can be expanded to other flow types
%                   'counter' for counter flow
%                   'parallel' for parallel
%   HeatExUCond     Condensor Overall heat transfer coefficient in W/(K m2)
%   HeatExACond     Condensor Heat transfer area in m2


% Outputs include two structures: ORC_Power and ORC_Temperature
% ORC_Power fields include:
%   outMech         Mechanical power output of expander in W
%   outHeat         Heat dissipation rate of condensor in W
%   inMech          Mechanical power input of pump in W
%   inHeat          Heat absorption rate of evaporator in W
% ORC_Temperature fields include:
%   source_in       Source fluid temperature at inlet of source heat exchanger in K
%   source_out      Source fluid outlet temperature in K
%   sink_in         Sink fluid temperature at inlet of sink heat exchanger in K
%   sink_out        Sink fluid outlet temperature in K
%   working_init    Working fluid initial temperature at inlet of source heat exchanger in K
%   expander_in     Working fluid temperature at inlet of expander in K
%   expander_out    Working fluid temperature at outlet of expander in K
%   pump_in         Working fluid temperature at inlet of pump in K
%   pump_out        Working fluid temperature at outlet of pump in K

% if no inputs are supplied then use the following values:
if nargin == 0
% Fluid name strings
    workingFluid = 'R245fa'; %'SES36';
    sourceFluid = 'Water';
    sinkFluid = 'Water';

% Mass flow rates of fluid in kg/s
    m_dot_working = 2.149;
    m_dot_source = 6 * 0.9619;%7.28;  % 12 kg/s is just shy of 200 gal/min of water
    m_dot_sink = 9 * 0.99997;%7.53;    % 12 kg/s is just shy of 200 gal/min of water

% Temperautres of source, sink, and working fluid in K
    T_source = 95 + 273.15;%(194.75+459.67)*5/9;%352.6;
    T_sink = 5 + 273.15;%(52.678+459.67)*5/9;%284.6;
    T_working_init = 295.6; %365.2; %12 + 273; %30 + 273;

% one atmosphere of pressure in pa for the pressure of the source and sink
    atm = 101325;

% High and low pressure of working fluid in Pa
    pressure_h = 6.0e5;%atm+5.8e5;%1.0e6;   %8.8e5; %
    pressure_l = 1.4e5; %atm+.55e5;%1.7e5; %0.6e5; %
    pressure_source = atm;
    pressure_sink = atm;

% pump and expander isentropic efficiencies
    eff_pump = 0.7;
    eff_expander = 0.78;
    
% heat exchanger parameters for evaporator
    heatExTypeVap = 'counter';
    heatExAVap = 66.5;%33.26; %9.3; %12; %268.5
    heatExUVap = 1500;
    
% heat exchanger parameters for condensor
    heatExTypeCond = 'counter';
    heatExACond = 66.5;%28.3; %18.6; %24; %268.5
    heatExUCond = 1400;
    
    plotflag = 1;
% % heat exchanger parameters -- Change to input
%     heatExType = 'counter';
%     heatExA = 5;
%     heatExU = 200;
else
    plotflag = 0;
end

% Initial Enthalpies of Source, Sink, and working fluid
H_source_in = CoolProp.PropsSI('H', 'P', pressure_source, 'T', T_source, sourceFluid);
H_sink_in = CoolProp.PropsSI('H', 'P', pressure_sink, 'T', T_sink, sinkFluid);
H_working_init = CoolProp.PropsSI('H', 'P', pressure_h, 'T', T_working_init, workingFluid);
% H_working_init = 233789.4;


% Initialize output structures to NaN
ORC_Power = struct( 'outMech', NaN, ...
                    'outHeat', NaN, ...
                    'inMech', NaN, ...
                    'inHeat', NaN);

ORC_Temperature = struct(   'source_in', NaN, ...
                            'source_out', NaN, ...
                            'sink_in', NaN, ...
                            'sink_out', NaN, ...
                            'working_init', NaN, ...
                            'expander_in', NaN, ...
                            'expander_out', NaN, ...
                            'pump_in', NaN, ...
                            'pump_out', NaN);

ORC_Temperature.source_in = T_source;
ORC_Temperature.sink_in = T_sink;
ORC_Temperature.working_init = T_working_init;
S_working_init = CoolProp.PropsSI('S', 'P', pressure_h, 'H', H_working_init, workingFluid);

% disp('working fluid phase before boiler:');
% disp(CoolProp.PhaseSI('P', pressure_h, 'T', T_working_init, workingFluid));

% Calculate the heat input at the source and the temperature changes 
[ ORC_Power.inHeat, H_source_out, H_expander_in ] = HeatExNTU_3(...
    H_source_in, m_dot_source, pressure_source, sourceFluid, H_working_init, m_dot_working, ...
    pressure_h, workingFluid, heatExTypeVap, heatExUVap, heatExAVap, 'H' );
% [ ORC_Power.inHeat, ORC_Temperature.source_out, ORC_Temperature.expander_in ] = HeatExNTU_2(...
%     T_source, m_dot_source, pressure_source, sourceFluid, T_working_init, m_dot_working, ...
%     pressure_h, workingFluid, heatExTypeVap, heatExUVap, heatExAVap);
% [ ORC_Power.inHeat, ORC_Temperature.expander_in, ORC_Temperature.source_out ] = HeatExNTU(...
%     T_source, m_dot_source, Cp_source, T_working_init, m_dot_working, ...
%     Cp_working, heatExType, heatExU, heatExA )

% % Used with HeatExNTU_3 to go from H to T
ORC_Temperature.source_out = CoolProp.PropsSI('T', 'P', pressure_source, 'H', H_source_out, sourceFluid);
ORC_Temperature.expander_in = CoolProp.PropsSI('T', 'P', pressure_h, 'H', H_expander_in, workingFluid);
S_expander_in = CoolProp.PropsSI('S', 'P', pressure_h, 'H', H_expander_in, workingFluid);

% power into the system is arbitrarily defined as negative
ORC_Power.inHeat = -ORC_Power.inHeat;

% disp('working fluid phase before expander:');
% disp(CoolProp.PhaseSI('P', pressure_h, 'H', H_expander_in, workingFluid));

% Calculate power output of isentropic change (expander)
[ORC_Power.outMech, H_expander_out] = IsentropicPower(H_expander_in, ...
    pressure_h, pressure_l, m_dot_working, workingFluid, eff_expander, 'H');
% [ORC_Power.outMech, ORC_Temperature.expander_out] = IsentropicPower(ORC_Temperature.expander_in, ...
%     pressure_h, pressure_l, m_dot_working, workingFluid, eff_expander, 'T');

% % Used with HeatExNTU_3 to go from H to T
ORC_Temperature.expander_out = CoolProp.PropsSI('T', 'P', pressure_l, 'H', H_expander_out, workingFluid);
S_expander_out = CoolProp.PropsSI('S', 'P', pressure_l, 'H', H_expander_out, workingFluid);

% disp('working fluid phase before condensor:');
% disp(CoolProp.PhaseSI('P', pressure_l, 'H', H_expander_out, workingFluid));

% Calculate the heat output at the sink and the temperature changes
[ ORC_Power.outHeat, H_pump_in, H_sink_out ] = HeatExNTU_3(...
    H_expander_out, m_dot_working, pressure_l, workingFluid, H_sink_in, m_dot_sink, ...
    pressure_sink, sinkFluid, heatExTypeCond, heatExUCond, heatExACond, 'H' );
% [ ORC_Power.outHeat, ORC_Temperature.pump_in, ORC_Temperature.sink_out ] = HeatExNTU_2(...
%     ORC_Temperature.expander_out, m_dot_working, pressure_l, workingFluid, ORC_Temperature.sink_in, m_dot_sink, ...
%     pressure_sink, sinkFluid, heatExTypeCond, heatExUCond, heatExACond);
% [ ORC_Power.outHeat, ORC_Temperature.pump_in, ORC_Temperature.sink_out ] = HeatExNTU(...
%     ORC_Temperature.expander_out, m_dot_working, Cp_working, T_sink, m_dot_sink, ...
%     Cp_sink, heatExType, heatExU, heatExA )

% % Used with HeatExNTU_3 to go from H to T
ORC_Temperature.pump_in = CoolProp.PropsSI('T', 'P', pressure_l, 'H', H_pump_in, workingFluid);
ORC_Temperature.sink_out = CoolProp.PropsSI('T', 'P', pressure_sink, 'H', H_sink_out, sinkFluid);
S_pump_in = CoolProp.PropsSI('S', 'P', pressure_l, 'H', H_pump_in, workingFluid);

% disp('working fluid phase before pump:');
% disp(CoolProp.PhaseSI('P', pressure_l, 'H', H_pump_in, workingFluid));

% Calculate power output of isentropic change (pump)
% ORC_Power.inMech = IsentropicPower(ORC_Temperature.pump_in, ORC_Temperature.pump_out, ...
[ORC_Power.inMech, H_pump_out ] = IsentropicPower(H_pump_in, ...
    pressure_l, pressure_h, m_dot_working, workingFluid, eff_pump, 'H');
% [ORC_Power.inMech, ORC_Temperature.pump_out ] = IsentropicPower(ORC_Temperature.pump_in, ...
%     pressure_l, pressure_h, m_dot_working, workingFluid, eff_pump, 'T');
% ORC_Power.inMech = abs(ORC_Power.inMech) %Return all Power values as positive for now

% % Used to go from H to T
ORC_Temperature.pump_out= CoolProp.PropsSI('T', 'P', pressure_h, 'H', H_pump_out, workingFluid);
S_pump_out = CoolProp.PropsSI('S', 'P', pressure_h, 'H', H_pump_out, workingFluid);

% disp('working fluid phase after pump:');
% disp(CoolProp.PhaseSI('P', pressure_h, 'H', H_pump_out, workingFluid));

if plotflag

    [pt_handle, hp_handle, st_handle] = R245faPTcurve(workingFluid);
    
    figure(pt_handle)
    hold on
    plot(T_working_init,pressure_h,'or')
    plot(ORC_Temperature.expander_in,pressure_h,'xm')
    plot(ORC_Temperature.expander_out,pressure_l,'xc')
    plot(ORC_Temperature.pump_in,pressure_l,'xg')
    plot(ORC_Temperature.pump_out,pressure_h,'xb')
    hold off
    
    figure(hp_handle)
    hold on
    plot(H_working_init,pressure_h,'or')
    plot(H_expander_in,pressure_h,'xm')
    plot(H_expander_out,pressure_l,'xc')
    plot(H_pump_in,pressure_l,'xg')
    plot(H_pump_out,pressure_h,'xb')
    hold off
    
    figure(st_handle)
    hold on
    plot(S_working_init,T_working_init,'or')
    plot(S_expander_in,ORC_Temperature.expander_in,'xm')
    plot(S_expander_out,ORC_Temperature.expander_out,'xc')
    plot(S_pump_in,ORC_Temperature.pump_in,'xg')
    plot(S_pump_out,ORC_Temperature.pump_out,'xb')
    hold off
    
end


end

