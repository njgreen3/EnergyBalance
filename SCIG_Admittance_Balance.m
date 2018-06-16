function [P_loss, Q_loss, Vph, f, P_mech ] = SCIG_Admittance_Balance( P_load, Q_load, n_mech, f_rated, Vline_rated, poles, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx )
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
    P_load = 4e6;
    Q_load = 0;
    f_rated = 60;
    poles = 6;
    n_mech = 1.01*f_rated*120/poles;
    Vline_rated = 6600;
    K_b = 0.5;
    K_w = 0.05;
    R1 = 0.044;
    L1 = 0.002;
    R2 = 0.078;
    L2 = 0.00029;
    Rm = 12000;
    Lm = 0.12;
    Rx = 0;
    Cx = 5.50e-05;
    
end

n_synch1d = linspace(0.99*f_rated*120/poles,1.01*n_mech,1024);
Vline1d = linspace(Vline_rated*0.8,Vline_rated*1.5,1024);

n_synch = repmat(n_synch1d', size(Vline1d));
Vline = repmat(Vline1d, size(n_synch1d'));
Vph = Vline/sqrt(3);

% calcluate slip. Should be negative when generating power
slip = (n_synch - n_mech)./n_synch;

% The electrical angular velocity w = 2pi*f = 2pi*P*n/120
w_e = 2*pi*poles.*n_synch/120;
f = w_e/2/pi;

% The mechanical angular velocity
w_mech = 2*pi*n_mech/60;

% calculate per-unit reactive impedance from inductance values (X = w*L)
% and capacitance value (X = 1/(w*C))
X1 = w_e.*L1;
X2 = w_e.*L2;
Xm = w_e.*Lm;
Xx = 1./(Cx*w_e);

% calculate friction and windage losses
P_bearing = K_b .* w_mech;
P_windage = K_w .* w_mech.^2;

% calculate stator, rotor, core, and excitation per-unit impedances
Z_stator = R1 + 1i*X1;   % R1 and X1 in series
Z_core = 1./(1./Rm + 1./(1i*Xm));   % Rm and Xm in parallel
Z_rotor = R2./slip + 1i*X2;   % R2/slip and X2 in series
Z_machine = Z_stator + 1./(1./Z_core + 1./Z_rotor);
Z_excite = Rx + Xx/1i;   %Rx and Xx in series
% Z_load = LoadImpedance(P_load, Q_load, Vline);  
Z_load = LoadImpedance(P_load, Q_load, Vline_rated);


Ymachine = 1./Z_machine;
Yexcite = 1./Z_excite;
Yload = 1./Z_load;

Smachine = Vline.^2.*conj(Ymachine);
Sexcite = Vline.^2.*conj(Yexcite);
Sload = Vline.^2.*conj(Yload);

if nargin == 0
    figure(1)
    surfc(f,Vph,real(Ymachine+Yexcite+Yload),'EdgeColor','flat')
    title('Frequency vs. Phase Voltage vs. Total Conductance')
    xlabel('F (Hz)')
    ylabel('V (V)')
    zlabel('G (siemens)')

    figure(2)
    surfc(f,Vph,imag(Ymachine+Yexcite+Yload),'EdgeColor','flat')
    title('Frequency vs. Phase Voltage vs. Total Susceptance')
    xlabel('F (Hz)')
    ylabel('V (V)')
    zlabel('B (siemens)')

    figure(3)
    surfc(f,Vph,abs(Ymachine+Yexcite+Yload),'EdgeColor','flat')
    title('Frequency vs. Phase Voltage vs. Total Admittance')
    xlabel('F (Hz)')
    ylabel('V (V)')
    zlabel('Y (siemens)')
    
    figure(4)
    surfc(f,Vph,real(Smachine+Sexcite+Sload),'EdgeColor','flat')
    title('Frequency vs. Phase Voltage vs. Total Conductance')
    xlabel('F (Hz)')
    ylabel('V (V)')
    zlabel('P (W)')

    figure(5)
    surfc(f,Vph,imag(Smachine+Sexcite+Sload),'EdgeColor','flat')
    title('Frequency vs. Phase Voltage vs. Total Susceptance')
    xlabel('F (Hz)')
    ylabel('V (V)')
    zlabel('Q (VAR)')

    figure(6)
    surfc(f,Vph,abs(Smachine+Sexcite+Sload),'EdgeColor','flat')
    title('Frequency vs. Phase Voltage vs. Total Admittance')
    xlabel('F (Hz)')
    ylabel('V (V)')
    zlabel('S (VA)')
    
    figure(7)
    plot(f(:,1),real(Ymachine(:,1)),f(:,1),real(Yexcite(:,1)),f(:,1),real(Yload*ones(size(f(:,1)))),f(:,1),real(Ymachine(:,1)+Yexcite(:,1)+Yload))
    title('Frequency vs. Conductance')
    xlabel('F (Hz)')
    ylabel('G (siemens)')
    legend('Machine','External Excitation','Load','Total','Location','best')

    figure(8)
    plot(f(:,1),imag(Ymachine(:,1)),f(:,1),imag(Yexcite(:,1)),f(:,1),imag(Yload*ones(size(f(:,1)))),f(:,1),imag(Ymachine(:,1)+Yexcite(:,1)+Yload))
    title('Frequency vs. Total Susceptance')
    xlabel('F (Hz)')
    ylabel('B (siemens)')
    legend('Machine','External Excitation','Load','Total','Location','best')

    figure(9)
    plot(f(:,1),abs(Ymachine(:,1)),f(:,1),abs(Yexcite(:,1)),f(:,1),abs(Yload*ones(size(f(:,1)))),f(:,1),abs(Ymachine(:,1)+Yexcite(:,1)+Yload))
    title('Frequency vs. Total Admittance')
    xlabel('F (Hz)')
    ylabel('Y (siemens)')
    legend('Machine','External Excitation','Load','Total','Location','best')

end

% find the index where sum of admittances about 0
errY = 0.005;
% errcount = 0;
% fvIDX = find(abs(Ymachine+Yexcite+Yload)<errY);
% while length(fvIDX) ~= 1
%     if isempty(fvIDX)
%         errY = errY*1.2;
%     else
%         errY = errY*0.8;
%     end
%     errcount = errcount + 1;
%     if errcount >= 1e3;
%         error('Could not balance admittances for rated values and parameters.')
%     end
%     fvIDX = find(abs(Ymachine+Yexcite+Yload)<errY);
% end
% disp([errcount, errY])

Ytotal = abs(Ymachine+Yexcite+Yload);
[minY, fvIDX] = min(Ytotal(:));
f = f(fvIDX);
Vph = Vph(fvIDX);
if minY > errY
    P_loss.mech = NaN;
    P_loss.copper = NaN;
    P_loss.core = NaN;
    P_loss.excitation = NaN;
    P_loss.total = NaN;
    Q_loss.rotor = NaN;
    Q_loss.stator = NaN;
    Q_loss.core = NaN;
    Q_loss.excitation = NaN;
    Q_loss.total = NaN;
    P_mech = NaN;
    if nargin == 0
        warning(['Could not balance admittances for rated values and parameters. Minimun Admittance of ' num2str(minY) ' siemens found at f = ' num2str(f) ' Hz and V = ' num2str(Vph) ' V.'])
    end
    Vph = NaN;
    f = NaN;
    return
end

Z_machine = Z_machine(fvIDX);
Z_excite = Z_excite(fvIDX);
Z_load = Z_load;%(fvIDX);
Z_stator = Z_stator(fvIDX);
Z_rotor = Z_rotor(fvIDX);
Z_core = Z_core(fvIDX);
X1 = X1(fvIDX);
X2 = X2(fvIDX);
Xm = Xm(fvIDX);
Xx = Xx(fvIDX);
 


% Calculate stator current 
I1 = Vph/Z_machine;

% Calculate excitation current
Ix = Vph./Z_excite;

E = Vph - I1*Z_stator;

% Calculate magnetizing current
Im = E./Z_core;

% Calculate rotor current
I2 = E./Z_rotor;

% Calculate load current
Iload = Vph./Z_load;

 
% Calculate per-unit active and reactive power losses in the asynch machine
Q_loss.rotor = 3*abs(I2).^2.*X2;
P_RCL = 3*abs(I2).^2.*R2;

Q_loss.stator = 3*abs(I1).^2.*X1;
P_SCL = 3*abs(I1).^2.*R1;

Q_loss.core = 3*abs(E).^2./Xm;
P_core = 3*abs(E).^2./Rm;

if any(Ix)
%     Sometimes when there is an infinite impedance matlab will calculate
%     NaN even though there is 0 current. This if/else ensures the power
%     calculations result in 0 when current is 0.
    Q_loss.excitation = -3*abs(Ix).^2.*Xx;
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

P_load = 3*abs(Vph)^2/conj(Z_load);
P_mech = P_load + P_loss.total;


