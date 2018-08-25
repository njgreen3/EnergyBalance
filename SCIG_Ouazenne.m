function [P_loss, Q_loss, Vph, f_e, P_mech, P_load_out ] = SCIG_Ouazenne( P_load_in, Q_load, n_mech, f_rated, Vline_rated, poles, K_b, K_w, R1, L1, R2, L2, magCurve, Rx, Cx )
%SCIG_ENERGY_BALANCE Summary of this function goes here 
%   This function calcluates the active power (P_loss) losses of a 3-phase 
%   squirrel cage induction generator (SCIG). Also calclulated is the
%   reactive power loss (Q_loss). Induction machines always consume Q
%   (negative Q_out) and the excitation capacitance generates Q (positive
%   Q_out).

%   Circuit diagram of an Asynchronous machine. 
%      Iin>   Ix         I1>        Im        I2>
%              v      R1     jw*X1   v     R2     jw*X2
%     0--------|-----RRRR----XXXX----|----RRRR----XXXX----|
%              |                    E|                    |
%     +        R  Rx               |-|-|                  |
%              R                   |   |                  R
%     Vph      |               Rm  R   X   jw*Xm          R  R2(1-slip)/slip
%              |                   R   X                  R
%     -       CCC  1/jw*Cx         |   |                  R
%              |                   |-|-|                  |
%              |                     |                    |  
%     0--------|---------------------|--------------------|
  
% Inputs:
%   P_load      Three phase active load power in W
%   Q_load      Three phase reactive load power in VAR
%   n_mech      Angular velocity of rotor in rpm
%   f_rated     Rated synchronous frequency in Hz
%   Vline_rated Rated line voltage in V
%   poles       Number of poles
%   K_b         Bearing coefficient 
%   K_w         Windage coefficient
%   R1          Stator resistance in ohm
%   L1          Stator inductance in H
%   R2          Rotor resistance in ohm
%   L2          Rotor inductance in H
%   Rm          Magnetizing resistance of the core in ohm
%   Lm          Magnetizing inductance of the core in L
%   Rx          Equivalent Series Resistance (ESR) of external excitation
%   Cx          External Excitation capacitance in F
%
% Outputs:
%   P_loss
%   Q_loss
%   Vph
%   P_mech      Mech power to rotor (P_mech = w_mech*Torque)

if nargin == 0
    P_load_in = 600;%7.46e3;
    Q_load = 0;
    f_rated = 60;
    poles = 2;
%     n_mech = 1.04*f_rated*120/poles;
    n_mech = 3000;
    Vline_rated = 220;
    K_b = 0.4;
    K_w = 0.001;
    R1 = 0.148;
    L1 = 0.423/(2*pi*f_rated);
    R2 = 0.144;
    L2 = 0.252/(2*pi*f_rated);
    [magCurve.X, magCurve.Vg_f] = imageinterpOuazenne;
    Rm = inf;
    Rx = 0;
    Cx = 240e-06;
    
end

n_synchRated = f_rated*120/poles;
b = n_mech/n_synchRated;

% The mechanical angular velocity
w_mech = 2*pi*n_mech/60;

% calculate reactive impedance from inductance values (X = w*L) and 
% capacitance value (X = 1/(w*C)) at rated frequencies
X1 = 2*pi*f_rated.*L1;
X2 = 2*pi*f_rated.*L2;
Xx = 1./(Cx*2*pi*f_rated);

% calculate friction and windage losses
P_bearing = K_b .* w_mech;
P_windage = K_w .* w_mech.^2;

% calculated the impedance of the load based off power consumed at rated
% line to line voltage
Z_load = LoadImpedance(P_load_in, Q_load, Vline_rated);
R_load = real(Z_load);

% Calcluate coefficients used to calculated the operating frequency
% relative to the rated frequency. These coefficients can be found in the
% appendix of Analysis of the Isolated Inducution Generator by Ouazenne and
% McPherson.  Note in the text Q2 and Q3 both have a typo which has been
% corrected here.
if R_load < inf
    Rtemp = R1 + R_load;
    
    Q0 = -b*R2*(Rtemp/R_load)^2;
    Q1 = R2*(Rtemp/R_load)^2 + Rtemp*(R2/R_load)^2 + b^2*Rtemp*(X2/R_load)^2;
    Q2 = -2*b*Rtemp*(X2/R_load)^2 - b*R2*((R1/Xx)^2 + (X1/R_load)^2 - 2*(X1/Xx));
    Q3 = R2*((X1/R_load)^2 + (R1/Xx)^2 - 2*(X1/Xx)) + Rtemp*(X2/R_load)^2 + R1*(R2/Xx)^2 + b^2*R1*(X2/Xx)^2;
    Q5 = R2*(X1/Xx)^2 + R1*(X2/Xx)^2;
    Q4 = -b*(Q5 + R1*(X2/Xx)^2);

else
    Q0 = -b*R2;
    Q1 = R2;
    Q2 = - b*R2*((R1/Xx)^2 - 2*(X1/Xx));
    Q3 = R2*((R1/Xx)^2 - 2*(X1/Xx)) + R1*(R2/Xx)^2 + b^2*R1*(X2/Xx);
    Q5 = R2*(X1/Xx)^2 + R1*(X2/Xx)^2;
    Q4 = -b*(Q5 + R1*(X2/Xx)^2);
end

% Calculate all the roots based from the coeffecients
aroots = roots([Q5 Q4 Q3 Q2 Q1 Q0]);

% Ouazzene and McPherson state there should only be one real root
a = aroots(~imag(aroots));

% Calculate the series equivalent of load and ext. excitation capacitance
RL = R_load*Xx^2/(a*(a^2*R_load^2 + Xx^2));
XL = R_load^2*Xx/(a^2*R_load^2 + Xx^2);

% calculate the operating frquency
f_e = a*f_rated;

% Calculate impedance of the core base of Ouazenne and McPherson
if R_load < inf
    Xm = ((RL + R1/a)*(X2^2 + (R2/(a-b))^2))/((R2/(a-b))*(X1-XL) - X2*(RL + R1/a));
else
    Xm = ((R1/a)*(X2^2 + (R2/(a-b))^2))/((R2/(a-b))*(X1 - Xx/a^2) - R1*X2/a);
end

% from Xm and magnatization curve calculate internal voltage magnitude
E_f = interp1(magCurve.X,magCurve.Vg_f,Xm,'pchip');
% calcluate phase voltage magnitude
if R_load < inf
    Vph = a*E_f/sqrt((1-a^2*(X1/Xx))^2 + a^2*((X1/R_load)^2 + (R1/Xx)^2) + (1 + R1/R_load)^2 - 1);
else
    Vph = a*E_f/sqrt((1-a^2*(X1/Xx))^2 + (a*R1/Xx)^2);
end

% calculate slip from electrical frequency and rotor speed relative to
% rated values
slip = (a-b)/a;

% % calculate stator, rotor, core, and excitation per-unit impedances
% Z_stator = R1/a + 1i*X1;   % R1 and X1 in series
% % Z_core = 1./(1./Rm + 1./(1i*Xm));   % Rm and Xm in parallel
% Z_core = 1i*Xm;
% Z_rotor = R2./(a*slip) + 1i*X2;   % R2/slip and X2 in series
% Z_machine = Z_stator + 1./(1./Z_core + 1./Z_rotor);
% Z_excite = Rx/a + Xx/(1i);   %Rx and Xx in series
% Z_lxs = Z_stator + 1/(1/Z_load + 1/Z_excite);

% Ymachine = 1./Z_machine
% Yexcite = 1./Z_excite
% Yload = 1./(R_load/a)
% 
% Y_rotor = 1/Z_rotor
% Y_core = 1/Z_core
% Y_lxs = 1/Z_lxs

% recalculate reactive impedances for operating frequency
X1 = 2*pi*f_e.*L1;
X2 = 2*pi*f_e.*L2;
Xx = 1./(Cx*2*pi*f_e);

% calculate stator, rotor, core, and excitation per-unit impedances
Z_stator = R1 + 1i*X1;   % R1 and X1 in series
% Z_core = 1./(1./Rm + 1./(1i*Xm));   % Rm and Xm in parallel
Z_core = 1i*Xm;
Z_rotor = R2./slip + 1i*X2;   % R2/slip and X2 in series
Z_machine = Z_stator + 1./(1./Z_core + 1./Z_rotor);
Z_excite = Rx + Xx/(1i);   %Rx and Xx in series
% Z_lxs = Z_stator + 1/(1/Z_load + 1/Z_excite);

% Ymachine = 1./Z_machine
% Yexcite = 1./Z_excite
% Yload = 1./Z_load
% 
% Y_rotor = 1/Z_rotor
% Y_core = 1/Z_core
% Y_lxs = 1/Z_lxs

% Calculate stator current 
I1 = Vph/Z_machine;

% Calculate excitation current
Ix = Vph./Z_excite;

% Calculate internal voltage
E = Vph - I1*Z_stator;

% Calculate magnetizing current
Im = E./Z_core;

% Calculate rotor current
I2 = E./Z_rotor;

% Calculate load current
Iload = Vph./Z_load;

 
% Calculate active and reactive power losses in the asynch machine
Q_loss.rotor = 3*abs(I2).^2.*X2;
P_RCL = 3*abs(I2).^2.*R2;

Q_loss.stator = 3*abs(I1).^2.*X1;
P_SCL = 3*abs(I1).^2.*R1;

Q_loss.core = 3*abs(E).^2./(Xm);
P_core =0;% 3*abs(E).^2./Rm;

if any(Ix)
%     Sometimes when there is an infinite impedance matlab will calculate
%     NaN even though there is 0 current. This if/else ensures the power
%     calculations result in 0 when current is 0.
    Q_loss.excitation = -3*abs(Ix).^2.*(Xx);
    P_excitation = 3*abs(Ix).^2.*Rx;
else
    Q_loss.excitation = 0;
    P_excitation = 0;
end

 
% Calculate total real and reactive powers. Reactive power is only
% consumed not generated therefore Q_out should be negative.
Q_loss.total = (Q_loss.core + Q_loss.rotor + Q_loss.stator);% + Q_loss.excitation;
% P_out = P_conv -(P_excitation + P_RCL + P_SCL + P_core);

% Output real losses in stucture
P_loss.mech = P_bearing + P_windage;
P_loss.copper = P_RCL + P_SCL;
P_loss.core = P_core;
P_loss.excitation = P_excitation;
P_loss.total = P_bearing + P_windage + P_RCL + P_SCL + P_core+ P_excitation;

P_load_out = 3*abs(Vph)^2/conj(Z_load);
P_mech = P_load_out + P_loss.total;


