function [ P_out, Q_out, P_loss, Q_loss] = SCIG_Energy_Balance_excitation( Vll, n_mech, f_synch, K_b, K_w, R1, X1, R2, X2, Rm, Xm, Rx, Xx )
%SCIG_ENERGY_BALANCE Summary of this function goes here
%   This function calcluates the active power (P_out) produced by a 3-phase 
%   squirrel cage induction generator (SCIG). Also calclulated is the
%   reactive power (Q_out) generated. Induction machines always consume Q
%   (negative Q_out) and the excitation capacitance generates Q (positive
%   Q_out).

%   Circuit diagram of an Asynchronous machine. 
%      Iin>   Ix         I1>        Im        I2>
%              v      R1     j*X1   v     R2     j*X2
%     0--------|-----RRRR----XXXX----|----RRRR----XXXX----|
%              |                    E|                    |
%     +        R  Rx               |-|-|                  |
%              R                   |   |                  R
%     Vph      |               Rm  R   X   j*Xm           R  R2(1-slip)/slip
%              |                   R   X                  R
%     -       CCC  1/jw*Cx         |   |                  R
%              |                   |-|-|                  |
%              |                     |                    |  
%     0--------|---------------------|--------------------|
  
% Inputs:
%   Vll         Line to Line voltage at the terminals of the machine
%   n_mech      Angular velocity of rotor in rpm
%   f_synch     Synchronous angular velocity in Hz
%   poles       Number of poles 
%   K_b         Bearing coefficient 
%   K_w         Windage coefficient
%   R1          Stator resistance in ohm
%   X1          Stator inductance in ohm
%   R2          Rotor resistance in ohm
%   X2          Rotor inductance in ohm
%   Rm          Magnetizing resistance of the core in ohm
%   Xm          Magnetizing inductance of the core in ohm
%   Rx          Equivalent Series Resistance (ESR) of external excitation
%   Xx          External Excitation capacitance in ohm
%
% Outputs:
%   P_out       Three phase net active power flow in W
%   Q_out       Three phase net reactive power flow in VAR
%   P_loss      Three phase active power consumed by the machine in W
%   Q_loss      Three phase reactive power consumed by the machine in VAR



if nargin == 0
    Vll = 220;
    f_synch = 60;
    n_mech = 914.8;
    poles = 4;
    K_b = 0.0;
    K_w = 0.0;
    R1 = 0.018;
    X1 = 0.32;
    R2 = 0.03;
    X2 = 0.05;
    Rm = 1800;
    Xm = 18;
    Rx = 0.003;
    Xx = 12;
    
%     n_synch1d = 1100:.5:n_mech-1;
%     Vline1d = 1e3:1e1:8e3;
%     n_synch = repmat(n_synch1d', size(Vline1d));
%     Vline = repmat(Vline1d, size(n_synch1d'));
%     Vph = Vline/sqrt(3);
end

% calculated the magnitude of the phase voltage. Phase angles will be
% relative to Vph
Vph = Vll/sqrt(3);

% Calculate synchronous and mechancial speeds in various units used later
% w_synch = 2*pi*f_synch;
n_synch = f_synch*120/(2*poles);
w_mech = 2*pi*n_mech/60;

% calcluate slip. Should be negative when generating power
slip = (n_synch - n_mech)./n_synch;

% calculate friction and windage losses
P_bearing = K_b .* w_mech;
P_windage = K_w .* w_mech.^2;


% calculate stator, rotor, core, and excitation per-unit impedances
Z_stator = R1 + 1i*X1;   % R1 and X1 in series
Z_core = 1./(1./Rm + 1./(1i*Xm));   % Rm and Xm in parallel
Z_rotor = R2./slip + 1i*X2;   % R2/slip and X2 in series
Z_machine = Z_stator + 1./(1./Z_core + 1./Z_rotor);
Z_excite = Rx + Xx/1i;   %Rx and Xx in series
Z_total = 1./(1./Z_excite + 1./Z_machine);
% Z_load = LoadImpedance(P_load, Q_load, Vline);


% Calculate stator current with thevenen equiv impedence of machine
I1 = Vph./Z_machine;

% Calculate internal voltage by subtracting drop across stator from Vph
E = Vph - I1.*Z_stator;

% Calculate refered rotor current
I2 = E./Z_rotor;

% Calculate excitation current
Ix = Vph./Z_excite;

% Calculate total real and reactive powers. 
% Calculate active and reactive power losses in the asynch machine
Q_rotor = 3*abs(I2).^2.*X2;
P_rotor = 3*abs(I2).^2.*R2;

Q_stator = 3*abs(I1).^2.*X1;
P_stator = 3*abs(I1).^2.*R1;

Q_core = 3*abs(E).^2./Xm;
P_core = 3*abs(E).^2./Rm;

Q_excitation = -3*abs(Ix).^2.*Xx;
P_excitation = 3*abs(Ix).^2.*Rx;


% Total losses/consumption
P_loss= P_bearing + P_windage + P_rotor + P_stator + P_core + P_excitation;
Q_loss = Q_core + Q_rotor + Q_stator;

% Net output (Produced - consumed)
P_out = -real(3*Vph^2./Z_total);
Q_out = -imag(3*Vph^2./Z_total);
% P_out = P_conv - (P_excitation + P_rotor + P_stator + P_core);


end

