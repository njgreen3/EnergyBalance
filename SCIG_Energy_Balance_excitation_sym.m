function [P_loss, Q_loss, Vph, E ] = SCIG_Energy_Balance_excitation( P_mech, P_load, Q_load, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx )
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
%   P_mech      Available per-unit Mech power to rotor (P_mech = w_mech*Torque)
%   P_load      Three phase per-unit active power 
%   Q_load      Three phase per-unit reactive power 
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
%   P_out
%   Q_out
%   P_loss
%   Vph
%   E

if nargin == 0
    P_mech = 1;
%     P_load = 0;
    P_load = .9;
    Q_load = 0;
    w_mech = 1.04;
    w_synch = 1;
    K_b = 0.01;
    K_w = 0.01;
%     V_leads = 1;
    R1 = 0.005;
    L1 = 0.088;
    R2 = 0.009;
    L2 = 0.125;
    Rm = 1400;
    Lm = 5;
    Rx = 0;
    Cx = 5;
end

% calcluate slip. Should be negative when generating power
slip = (w_synch - w_mech)./w_synch;

% calculate per-unit reactive impedance from inductance values (X = w*L)
% and capacitance value (X = 1/(w*C))
X1 = w_synch.*L1;
X2 = w_synch.*L2;
Xm = w_synch.*Lm;
Xx = 1./(Cx.*w_synch);

% calculate friction and windage losses
P_bearing = K_b .* w_mech;
P_windage = K_w .* w_mech.^2;

% remaining power is converted to electrical
P_conv = P_mech - P_bearing - P_windage;

% Calculate rotor current based off of slip and mechanical power after friction
I2 = abs(sqrt((P_conv/R2 .* slip./(1-slip))));

% calculate stator, rotor, core, and excitation per-unit impedances
Z_stator = R1 + 1i*X1;   % R1 and X1 in series
Z_core = 1./(1./Rm + 1./(1i*Xm));   % Rm and Xm in parallel
Z_rotor = R2./slip + 1i*X2   % R2/slip and X2 in series
Z_excite = Rx + Xx/1i;   %Rx and Xx in series

angle(Z_rotor)*180/pi
cos(angle(Z_rotor))

S_load = P_load + 1i*Q_load;

syms x;
% vph_eqn = x == (I2.*Z_excite.*Z_core*x*conj(x)./conj(S_load))./(x.*conj(x)./conj(S_load).*(Z_stator + Z_core + Z_excite) + Z_stator.*Z_excite + Z_core.*Z_excite);
vph_eqn = 1 == (I2.*Z_core.*conj(x)./conj(S_load))./(x.*conj(x)./conj(S_load).*(Z_stator./Z_excite + Z_core./Z_excite + 1) + Z_stator + Z_core);

Vph = eval(solve(vph_eqn, x))
abs(Vph)
angle(Vph)*180/pi

% % Initially assume load is across a line voltage of 0 per-unit to start while loop
% V_load = -inf;
% % Initially assume phase voltage is 1 to start while loop
% Vph = inf;
% % accectable error between assumed load across load and calculated Vph
% err = 0.00000001;
% % initialize count variable to break if the voltages do not converge
% vCount = 0;
% vLimit = 1e3;
% 
% % while( any(any(abs(V_load - abs(Vph)*sqrt(3)) > err)))
% while( any(any(abs(abs(Vph)*sqrt(3) - V_load) > err)))
% % Each cycle sets the load voltage to the magnitude of the previously 
% % calucated phase voltage (correcting for line to phase) and recalculate
% % equivalent impedance and recalculate
% 
% %   Check if voltage is taking too long to converge
%     if vCount >= vLimit
%         error(['Voltage did not converge after ' num2str(vCount) ' iterations.'])
%     end
% %   Increase count of cycle
%     vCount = vCount + 1;
%     
% %   Set V_load to magnitude of previously calulated value of Vph with line 
% %   to phase conversion
    V_load = abs(Vph);
%     
% %   Calculate equivalent load impedance for given line voltage    
    Z_load = LoadImpedance(P_load, Q_load, V_load);
% %   Calculate parallel combination of load impedance and excitaion impedance
    Z_load_par_excite = 1./(1./Z_load + 1./Z_excite);
% 
% %   Calculate internal and phase voltages of the circuit with angles 
% %   relative rotor current
    E = I2./(1./Z_core + 1./(Z_stator + Z_load_par_excite))
    abs(E)
    angle(E)*180/pi
%     Vph = E.*Z_load_par_excite./(Z_stator + Z_load_par_excite);
%     
% end
% Vph


% Calculate stator current with current divider
I1 = I2.*Z_core./(Z_core + Z_stator + Z_load_par_excite);

% Calculate excitation current
Ix = Vph./Z_excite;

% Calculate magnetizing current
Im = E./Z_core;

% Calculate load current
Iload = Vph./Z_load;

% if isinf(Xx)
%     I1 = 0;
%     Ix = 0;
%     Iload = 0;
% end

% Determine the phase offset between I2 and Vph in order to calculate the
% currents and voltages with phases relative to Vph
ang_offset = exp(-1i*angle(Vph));
if isnan(ang_offset)
    ang_offset = 1;
end

% Rotate current and voltage phases so they are all relatve to the phase
% voltage. Not really necessary here, but it may be helpful later.
Vph = Vph.*ang_offset;
E = E.*ang_offset;
Iload = Iload.*ang_offset;
Ix = Ix.*ang_offset;
I1 = I1.*ang_offset;
Im = Im.*ang_offset;
I2 = I2.*ang_offset;


% Calculate per-unit active and reactive power losses in the asynch machine
Q_loss.rotor = abs(I2).^2.*X2;
P_RCL = abs(I2).^2.*R2;

Q_loss.stator = abs(I1).^2.*X1;
P_SCL = abs(I1).^2.*R1;

Q_loss.core = abs(E).^2./Xm;
P_core = abs(E).^2./Rm;

if any(Ix)
%     Sometimes when there is an infinite impedance matlab will calculate
%     NaN even though there is 0 current. This if/else ensures the power
%     calculations result in 0 when current is 0.
    Q_loss.excitation = -abs(Ix).^2.*Xx;
    P_excitation = abs(Ix).^2.*Rx;
else
    Q_loss.excitation = 0;
    P_excitation = 0;
end
% if any(Iload)
% %     Sometimes when there is an infinite impedance matlab will calculate
% %     NaN even though there is 0 current. This if/else ensures the power
% %     calculations result in 0 when current is 0.
%     Q_out1 = 3*abs(Iload).^2.*imag(Z_load);
%     P_out1 = 3*abs(Iload).^2.*real(Z_load);
% else
%     Q_out1 = 0;
%     P_out1 = 0;
% end


% 3*abs(I2).^2*R2./slip ;
% 3*abs(I2).^2*R2./slip+P_SCL + P_core;

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

if P_mech < P_load
    warning('The load exceeds the available mechanical power.')
elseif P_mech < P_load + P_loss.total
    warning('The load and losses exceed the available mechanical power.')
end
end

