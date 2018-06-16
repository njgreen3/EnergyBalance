function [ P_out, Q_out, P_loss, Vph, E ] = SCIG_Energy_Balance( P_mech, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm )
%SCIG_ENERGY_BALANCE Summary of this function goes here
%   This function calcluates the active power (P_out) produced by a 3-phase 
%   squirrel cage induction generator (SCIG). Also calclulated is the
%   reactive power (Q_out) generated, but induction machines always consume
%   Q therefore Q_out should be negative.

%   Circuit diagram of an Asynchronous machine. 
%              I1>        Im        I2>
%           R1     w*X1    v     R2     w*X2
%     0----RRRR----XXXX----|----RRRR----XXXX----|
%                         E|                    |
%     +                  |-|-|                  |
%                        |   |                  R
%     Vph            Rm  R   X   w*Xm           R  R(1-slip)/slip
%                        R   X                  R
%     -                  |   |                  R
%                        |-|-|                  |
%                          |                    |  
%     0--------------------|--------------------|
  
% Inputs:
%   P_mech      Per-unit Mech power to rotor (P_mech = w_mech*Torque)
%   w_mech      Per-unit Angular velocity of rotor 
%   w_synch     Per-unit Synchronous angular velocity
%   K_b         Bearing coefficient 
%   K_w         Windage coefficient
%   V_leads     Line-Line per-unit voltage at the leads 
%   R1          Stator per-unit resistance
%   L1          Stator per-unit inductance
%   R2          Rotor per-unit resistance
%   L2          Rotor per-unit inductance
%   Rm          Magnetizing per-unit resistance of the core
%   Lm          Magnetizing per-unit inductance of the core

if nargin == 0
    P_mech = 1.0;
    w_mech = 1.05;
    w_synch = 1;
    K_b = .01;
    K_w = .01;
%     V_leads = 1;
    R1 = .005;
    L1 = .088;
    R2 = .009;
    L2 = .0125;
    Rm = 1400;
    Lm = 5;
end

% calcluate slip. Should be negative when generating power
slip = (w_synch - w_mech)/w_synch;

% calculate per-unit reactive impedance from inductance values (X = w*L)
X1 = w_synch*L1;
X2 = w_synch*L2;
Xm = w_synch*Lm;

% calculate friction and windage losses
P_bearing = K_b * w_mech;
P_windage = K_w * w_mech.^2;

% remaining power is converted to electrical
P_conv = P_mech - P_bearing - P_windage;

% calculate stator, rotor, core, and overall thevenin per-unit impedances
Z_stator = R1 + 1i*X1;   % R1 and X1 in series
Z_core = 1i*Rm*Xm/(Rm + 1i*Xm);   % Rm and Xm in parallel
Z_rotor = R2./slip + 1i*X2;   % R2/slip and X2 in series
Z_thev = Z_stator + Z_core.*Z_rotor./(Z_core + Z_rotor);  %Zcore and Zrotor in parallel and Zstator in series

% Calculate currents and voltages of the circuit with phases relative to I2
% from the diagram above all in per-unit 
I2 = abs(sqrt((P_conv/R2/3 .* slip./(1-slip))));
E = I2.*Z_rotor;
Im = E/Z_core;
I1 = Im + I2;
Vph = I1.*Z_thev;

% I2_ang = angle(I2)*180/pi;
% E_ang = angle(E)*180/pi;
% Im_ang = angle(Im)*180/pi;
% I1_ang = angle(I1)*180/pi;
% Vph_ang = angle(Vph)*180/pi;


% Determine the phase offset between I2 and V in order to calculate the
% currents and voltages with phases relative to Vph
ang_offset = exp(-1i*angle(Vph));

% Rotate current and voltage phases so they are all relatve to the phase
% voltage. Not really necessary here, but it may be helpful later.
Vph = Vph.*ang_offset;
I1 = I1.*ang_offset;
E = E.*ang_offset;
Im = Im.*ang_offset;
I2 = I2.*ang_offset;

% I2_ang = angle(I2)*180/pi;
% E_ang = angle(E)*180/pi;
% Im_ang = angle(Im)*180/pi;
% I1_ang = angle(I1)*180/pi;
% Vph_ang = angle(Vph)*180/pi;
% abs(E);

% Calculate per-unit active and reactive power losses in the asynch machine
Q_rotor = 3*abs(I2).^2*X2;
P_RCL = 3*abs(I2).^2*R2;

Q_core = 3*abs(E).^2/Xm;
P_core = 3*abs(E).^2/Rm;

Q_stator = 3*abs(I1).^2*X1;
P_SCL = 3*abs(I1).^2*R1;

% Calculate output real and reactive powers. Reactive power is only
% consumed not generated therefore Q_out should be negative.
Q_out = -(Q_core + Q_rotor + Q_stator);
P_out = P_conv -(P_RCL + P_SCL + P_core);

% Output real losses in stucture
P_loss.mech = P_bearing + P_windage;
P_loss.copper = P_RCL + P_SCL;
P_loss.core = P_core;


end

