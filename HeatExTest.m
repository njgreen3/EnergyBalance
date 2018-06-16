clear 
clc

% define some properties/variables for the heat exchanger itself
HeatExCoeff = 200;
HeatExArea = 5;
HeatExFlowtype = 'counter';

% Temerpature profile of heat source input in degC
T_hot_in = [90 90 91 90 89 90 91 92 93 94];

% Temerpature profile of heat sink input in degC
T_cool_in = [31 30 30 31 30 29 30 31 32 31];


% Specific heat values
% Water (liquid) properties
H2O_C = 4180; %J/(kg K)

% R227ea properties
R227_C = 820.6;  %J/(kg K)

% Flow values
% Heat source
V_dot_h = 200; %gal/min
% gallon to m3 conversion
gal_mcube = 264.172; %gal/m3
% min to sec conversion
min_sec = 60; %s/min
% convert volume flow rate of heat source from gal/min to m3/s
V_dot_conv = V_dot_h/gal_mcube/min_sec;
% Water density is roughly linear from 50 degC to 100 degC
% The following relation is mostly invalid outside this range
H2O_density = 1016 - (29*T_hot_in)/50; %kg/m3
% Calculate mass flow rate from density and volume flow rate
m_dot_h = V_dot_conv*H2O_density; %kg/s

%Heat Sink;
m_dot_c = 10; %kg/s

%Calculte the heat flow rate and the outlet temperatures
[Heat_rate, T_hot_out, T_cool_out] = HeatExNTU(T_hot_in, m_dot_h, H2O_C, ...
    T_cool_in, m_dot_c, R227_C, HeatExFlowtype, HeatExCoeff, HeatExArea)
%%
% hot and cool fluid names
fluid_h = 'water';
fluid_c = 'water';

% Hot and cool inlet temp in K
Thi = 90+273;
Tci = 330;

% Hot and cool mass flow rates in kg/s
m_dot_h = 12;   % just shy of 200 gal/min of water
m_dot_c = 1;

% Hot and cool pressures
p_h = 101325;   % about 1 atm at sea level
p_c = 5e5;

% define some properties/variables for the heat exchanger itself
U = 200;
A = 5;
HeatExType = 'counter';

[ P1, Tho1, Tco1 ] = HeatExNTU_2( Thi, m_dot_h, p_h, fluid_h,...
    Tci, m_dot_c, p_c, fluid_c, ...
    HeatExType, U, A );

[ P2a, Tho2a, Tco2a ] = HeatExNTU_2( Thi, m_dot_h, p_h, fluid_h,...
    Tci, m_dot_c, p_c, fluid_c, ...
    HeatExType, U, A/2 );

[ P2b, Tho2b, Tco2b ] = HeatExNTU_2( Tho2a, m_dot_h, p_h, fluid_h,...
    Tco2a, m_dot_c, p_c, fluid_c, ...
    HeatExType, U, A/2 );

[Tho1 Tco1; Tho2a Tco2a; Tho2b Tco2b;]
%%
clc
% Test for HeatExNTU_3
% hot and cool fluid names
fluid_h = 'water';
fluid_c = 'R245fa';

% Hot and cool inlet temp in K
Thi = 90+273;
Tci = 320;

% Hot and cool mass flow rates in kg/s
m_dot_h = 12;   % just shy of 200 gal/min of water
m_dot_c = 3;

% Hot and cool pressures
p_h = 101325;   % about 1 atm at sea level
p_c = 5e5;

Hhi = CoolProp.PropsSI('H','P',p_h,'T',Thi,fluid_h);
Hci = CoolProp.PropsSI('H','P',p_c,'T',Tci,fluid_c);


% define some properties/variables for the heat exchanger itself
U = 1500;
A = 10;
HeatExType = 'parallel';
warning('off','SWIG:RuntimeError')
disp('one')
[ P1, Hho1, Hco1 ] = HeatExNTU_3( Hhi, m_dot_h, p_h, fluid_h,...
    Hci, m_dot_c, p_c, fluid_c, ...
    HeatExType, U, A, 'H' );
disp('two a')
[ P2a, Hho2a, Hco2a ] = HeatExNTU_3( Hhi, m_dot_h, p_h, fluid_h,...
    Hci, m_dot_c, p_c, fluid_c, ...
    HeatExType, U, A/2, 'H' );
disp('two b')
[ P2b, Hho2b, Hco2b ] = HeatExNTU_3( Hho2a, m_dot_h, p_h, fluid_h,...
    Hco2a, m_dot_c, p_c, fluid_c, ...
    HeatExType, U, A/2, 'H' );
warning('on','SWIG:RuntimeError')

Tho1 = CoolProp.PropsSI('T','P',p_h,'H',Hho1,fluid_h);
Tco1 = CoolProp.PropsSI('T','P',p_c,'H',Hco1,fluid_c);
Tho2a = CoolProp.PropsSI('T','P',p_h,'H',Hho2a,fluid_h);
Tco2a = CoolProp.PropsSI('T','P',p_c,'H',Hco2a,fluid_c);
Tho2b = CoolProp.PropsSI('T','P',p_h,'H',Hho2b,fluid_h);
Tco2b = CoolProp.PropsSI('T','P',p_c,'H',Hco2b,fluid_c);

[Tho1 Tco1; Tho2a Tco2a; Tho2b Tco2b;]

%%
clc
% hot and cool fluid names
fluid_h = 'water';
fluid_c = 'R245fa';
% fluid_c = 'water';

% Hot and cool inlet temp in K
T_or_H_hi = (90:95)'+273;
T_or_H_ci = 330*ones(size(T_or_H_hi));

% Hot and cool mass flow rates in kg/s
m_dot_h = 12*ones(size(T_or_H_hi));   % just shy of 200 gal/min of water
m_dot_c = 2*ones(size(T_or_H_hi));

% Hot and cool pressures
p_h = 101325*ones(size(T_or_H_hi));   % about 1 atm at sea level
p_c = 4.5e5*ones(size(T_or_H_hi));

% define some properties/variables for the heat exchanger itself
U = 1500;
A = 10;
HeatExType = 'counter';
T_or_H = 'T';

[ P1, Tho1, Tco1 ] = HeatExNTU_3_low( T_or_H_hi, m_dot_h, p_h,  ...
    T_or_H_ci, m_dot_c, p_c, fluid_h, fluid_c, HeatExType, U, A, T_or_H )
% P1 = [timewithTh.Time P1];
% Tho1 = [timewithTh.Time Tho1];
% Tco1 = [timewithTh.Time Tco1];

% [ P2a, Tho2a, Tco2a ] = HeatExNTU_3_low( T_or_H_hi, m_dot_h, p_h, ...
%     T_or_H_ci, m_dot_c, p_c, fluid_h, fluid_c, HeatExType, U, A/2, T_or_H)
% 
% [ P2b, Tho2b, Tco2b ] = HeatExNTU_3_low( Tho2a, m_dot_h, p_h, ...
%     Tco2a, m_dot_c, p_c, fluid_h, fluid_c, HeatExType, U, A/2, T_or_H )

% [Tho1 Tho2a Tho2b Tco1 Tco2a Tco2b;]
% [P1, P2a+P2b]