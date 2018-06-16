function [P_loss, Q_loss, Vph, P_mech ] = SCIG_Energy_Balance_excitation( P_load, Q_load, n_mech, n_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx )
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
%   P_load      Three phase active load power in W
%   Q_load      Three phase reactive load power in VAR
%   n_mech      Angular velocity of rotor in rpm
%   n_synch     Synchronous angular velocity in rpm
%   pole        Number of poles
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
    P_load = 5e6;
    Q_load = 0;
    n_synch = 1200;
    n_mech = n_synch*1.04;
    poles = 6;
    K_b = 0.01;
    K_w = 0.01;
    R1 = 0.0436;
    L1 = 0.002;
    R2 = 0.0784;
    L2 = 0.000289;
    Rm = 12197;
    Lm = 0.1155;
    Rx = 0;
    Cx = 1e-4;
    
    n_synch1d = 1100:.5:n_mech-1;
    Vline1d = 1e3:1e1:8e3;
    n_synch = repmat(n_synch1d', size(Vline1d));
    Vline = repmat(Vline1d, size(n_synch1d'));
    Vph = Vline/sqrt(3);
end

% calcluate slip. Should be negative when generating power
slip = (n_synch - n_mech)./n_synch;

% The electrical angular velocity w = 2pi*f = 2pi*P*n/120
w_e = 2*pi*poles*n_synch/120;
f = w_e/2/pi;

% The mechanical angular velocity
w_mech = 2*pi*n_mech/60;

% calculate per-unit reactive impedance from inductance values (X = w*L)
% and capacitance value (X = 1/(w*C))
X1 = w_e.*L1;
X2 = w_e.*L2;
Xm = w_e.*Lm;
Xx = 1./(Cx.*w_e);

% calculate friction and windage losses
P_bearing = K_b .* w_mech;
P_windage = K_w .* w_mech.^2;

% calculate stator, rotor, core, and excitation per-unit impedances
Z_stator = R1 + 1i*X1;   % R1 and X1 in series
Z_core = 1./(1./Rm + 1./(1i*Xm));   % Rm and Xm in parallel
Z_rotor = R2./slip + 1i*X2;   % R2/slip and X2 in series
Z_machine = Z_stator + 1./(1./Z_core + 1./Z_rotor);
Z_excite = Rx + Xx/1i;   %Rx and Xx in series
Z_load = LoadImpedance(P_load, Q_load, Vline);

Ymachine = 1./Z_machine;
Yexcite = 1./Z_excite;
Yload = 1./Z_load;

fvIDX = find(abs(Ymachine+Yexcite+Yload)<.0004);
f(fvIDX)
Vph(fvIDX)
Z_machine(fvIDX)
Z_excite(fvIDX)
Z_load(fvIDX)
Vph(fvIDX).^2./conj(Z_machine(fvIDX))
Vph(fvIDX).^2./conj(Z_excite(fvIDX))
Vph(fvIDX).^2./conj(Z_load(fvIDX))

figure(1)
surfc(f,Vph,real(Ymachine+Yexcite+Yload),'EdgeColor','flat')

figure(2)
surfc(f,Vph,imag(Ymachine+Yexcite+Yload),'EdgeColor','flat')

figure(3)
surfc(f,Vph,abs(Ymachine+Yexcite+Yload),'EdgeColor','flat')


% Initially assume current differences is inf to start while loop
iDiff = inf;
% accectable error between assumed load across load and calculated Vph
err = 1e-9;
% initialize count variable to break if the currents do not converge
iCount = 0;
iLimit = 1e4;

% figure()
% plot(vCount, abs(Vph), '*')
% plot(vCount, abs(V_load/sqrt(3) - abs(Vph)), '*')
% hold on

while(any(any(abs(iDiff) > err)))
% Each cycle sets the load voltage to the magnitude of the previously 
% calucated phase voltage (correcting for line to phase) and recalculate
% equivalent impedance and recalculate

%     abs(V_load - abs(Vph))
    
%   Check if voltage is taking too long to converge
    if iCount >= iLimit
        error(['Current did not converge after ' num2str(iCount) ' iterations.'])
    end
%   Increase count of cycle
    iCount = iCount + 1;
    
%     plot(vCount, abs(Vph), '*')
%     plot(vCount, abs(V_load - abs(Vph)), '*')
    
%   Set V_load to magnitude of previously calulated value of Vph with line 
%   to phase conversion
    V_load = abs(Vph);
    
    
%   Calculate equivalent load impedance for given line voltage    
    Z_load = LoadImpedance(P_load, Q_load, V_load);
%   Calculate parallel combination of load impedance and excitaion impedance
    Z_load_par_excite = 1./(1./Z_load + 1./Z_excite);
%   Calculate internal and phase voltages of the circuit with angles 
%   relative rotor current
    E = I2./(1./Z_core + 1./(Z_stator + Z_load_par_excite));
    Vph = E.*Z_load_par_excite./(Z_stator + Z_load_par_excite);
%     abs(Vph)
    abs(V_load - abs(Vph));
    
end
% vCount
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
disp(angle(E)*180/pi);
Iload = Iload.*ang_offset;
Ix = Ix.*ang_offset;
I1 = I1.*ang_offset;
Im = Im.*ang_offset;
I2 = I2.*ang_offset;
disp(angle(I2)*180/pi);

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
%     Q_out1 = abs(Iload).^2.*imag(Z_load);
%     P_out1 = abs(Iload).^2.*real(Z_load);
% else
%     Q_out1 = 0;
%     P_out1 = 0;
% end


% abs(I2).^2*R2./slip ;
% abs(I2).^2*R2./slip+P_SCL + P_core;

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

