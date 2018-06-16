function [P_loss, Q_loss, P_mech, S_grid] = SCIG_Energy_Balance_excitation( P_load, Q_load, Vph, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx )
%SCIG_ENERGY_BALANCE Summary of this function goes here
%   This function calcluates the active power (P_out) produced by a 3-phase 
%   squirrel cage induction generator (SCIG). Also calclulated is the
%   reactive power (Q_out) generated. Induction machines always consume Q
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
%   P_load      Three phase per-unit active power 
%   Q_load      Three phase per-unit reactive power 
%   Vph         Terminal phase voltage in per-unit
%   w_mech      Per-unit Angular velocity of rotor 
%   w_synch     Per-unit Synchronous angular velocity
%   K_b         Bearing coefficient 
%   K_w         Windage coefficient
%   R1          Stator per-unit resistance
%   L1          Stator per-unit inductance
%   R2          Rotor per-unit resistance
%   L2          Rotor per-unit inductance
%   Rm          Magnetizing per-unit resistance of the core
%   Lm          Magnetizing per-unit inductance of the core
%   Rx          Per-unit Equivalent Series Resistance (ESR) of excitation
%   Cx          Per-unit excitation capacitance
%
% Outputs:
%   P_loss      Per-unit active power losses
%    mech           mechanical losses
%    rotor          copper losses in rotor
%    stator         copper losses in stator
%    core           iron losses in core
%    excite         ESR losses of external excitation brach
%    total          Sum of listed losses
%   Q_loss      Per-unit reactive power losses
%    rotor          Q consumed by rotor
%    stator         Q consumed by stator
%    core           Q consumed by core for excitation
%    excite     Q from external excitation brach (should be negative)
%    total          Sum of rotor, stator and core losses. This value
%                   excludes the Q provided by the external excitation
%   P_mech      Per-unit Mech power to shaft in order to achieve input
%                   conditions
%   S_grid      Remaining apparent power provided or consumed by the grid 
%                   or inverter in order to achieve input conditions.
%                   Negative values indicate the grid or inverter must
%                   consume or dissipate the excess in order to achieve 
%                   input conditions. Positive values indicate the grid or 
%                   inverter must provide the deficit in order to achieve
%                   input conditions.
%                   

if nargin == 0
%     P_mech = 1;
    P_load = .9;
    Q_load = 0;
    Vph = 1;
    w_mech = 1.01;
    w_synch = 1;
    K_b = 0.01;
    K_w = 0.01;
    R1 = 0.005;
    L1 = 0.088;
    R2 = 0.009;
    L2 = 0.0125;
    Rm = 1400;
    Lm = 5;
    Rx = .001;
    Cx = .3;
end

% calcluate slip. Should be negative when generating power
slip = (w_synch - w_mech)./w_synch;

% calculate per-unit reactive impedance from inductance values (X = w*L)
% and capacitance value (X = 1/(w*C))
X1 = w_synch.*L1;
X2 = w_synch.*L2;
Xm = w_synch.*Lm;
Xx = 1./(Cx.*w_synch);


% calculate stator, rotor, core, and excitation per-unit impedances
Z_stator = R1 + 1i*X1;   % R1 and X1 in series
Z_core = 1./(1./Rm + 1./(1i*Xm));   % Rm and Xm in parallel
Z_rotor = R2./slip + 1i*X2;   % R2/slip and X2 in series
Z_excite = Rx + Xx/1i;   %Rx and Xx in series
Z_load = LoadImpedance(P_load, Q_load, Vph);

I_load = Vph./Z_load;
I_excite = Vph./Z_excite;
I_stator = Vph./(Z_stator + (Z_core.^-1 + Z_rotor.^-1).^-1);

E = Vph - I_stator.*Z_stator;

I_core = E./Z_core;
I_rotor = E./Z_rotor;

P_conv = abs(abs(I_rotor).^2 .* R2.*(1-slip)./slip);

% calculate mechanical friction and windage losses
P_bearing = K_b .* w_mech;
P_windage = K_w .* w_mech.^2;
P_loss.mech = P_bearing + P_windage;

P_mech = P_conv + P_bearing + P_windage;

% Calculate per-unit active and reactive power losses in the asynch machine
Q_loss.rotor = abs(I_rotor).^2.*X2;
P_loss.rotor = abs(I_rotor).^2.*R2;

Q_loss.stator = abs(I_stator).^2.*X1;
P_loss.stator = abs(I_stator).^2.*R1;

Q_loss.core = abs(E).^2./Xm;
P_loss.core = abs(E).^2./Rm;

Q_loss.excite = -abs(I_excite).^2.*Xx;
P_loss.excite = abs(I_excite).^2.*Rx;

% Calculate total real and reactive powers. 
Q_loss.total = (Q_loss.core + Q_loss.rotor + Q_loss.stator);% + Q_loss.excite;
P_loss.total = P_loss.mech + P_loss.rotor + P_loss.stator + P_loss.core + P_loss.excite;

S_grid = Vph.*conj(I_load + I_excite + I_stator);


end

