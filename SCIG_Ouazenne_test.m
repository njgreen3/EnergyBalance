clear 
close all
clc 
%%
% This script tests the SCIG_Admittance_Balance function under a range of 
% inputs inculding external excitation capacitance, rotor speeds, and
% loads.

% 3 phase electrical load
P_load_in = 2e3:2e3:6e3;   %W
Q_load = 0;         %VAR
% rated values
f_rated = 60;       %Hz
Vline_rated = 220; %V
poles = 2;
% speeds
% n_mech = 1.01*f_rated*120/poles;
n_mech = 3000;
% friction coefficients
K_b = 0.4;
K_w = 0.001;
% stator values
R1 = 0.148;         %ohm
L1 = 0.423/2/pi/f_rated;     %H
% rotor values
R2 = 0.144;         %ohm
L2 = 0.252/2/pi/f_rated;       %H
% core values
[magCurve.X, magCurve.Vg_f] = imageinterp;
% Rm = 12000;         %ohm
% Lm = 0.12;          %H
% excitation values

Rx = 0;             %ohm
Cx_arr = linspace(1e-6,400e-6,500);   %F
% Cx = repmat(Cx,1,length(w_mech));

Vph_cx = zeros(size(Cx_arr));
f_cx = zeros(size(Cx_arr));
P_mech_cx = zeros(size(Cx_arr));
P_load_out_cx = zeros(size(Cx_arr));

% Test for different loads
for m = 1:length(P_load_in)
% Test varying capacitance
for n = 1:length(Cx_arr);
% test with varying w but constant P_mech
[P_loss_cx(m,n), Q_loss_cx(m,n), Vph_cx(m,n), f_cx(m,n), P_mech_cx(m,n), P_load_out_cx(m,n) ] = SCIG_Ouazenne( P_load_in(m), Q_load, n_mech, f_rated, Vline_rated, poles, K_b, K_w, R1, L1, R2, L2, magCurve, Rx, Cx_arr(n) );



end
end
unboundFlag = abs(Vph_cx) > 1e3 | Vph_cx < 0;
Vph_cx(unboundFlag) = NaN;
f_cx(unboundFlag) = NaN;
P_mech_cx(unboundFlag) = NaN;
P_load_out_cx(unboundFlag) = NaN;
% [Q_loss(unboundFlag).total] = [NaN];
% Q_loss(unboundFlag).excitation = NaN;

figure(1)
plot(   Cx_arr, P_mech_cx(1,:),Cx_arr, P_mech_cx(2,:),Cx_arr, P_mech_cx(3,:),...
        Cx_arr, P_load_out_cx(1,:),Cx_arr, P_load_out_cx(2,:),Cx_arr, P_load_out_cx(3,:))
title('Capacitance vs. Mech & Elec Power ')
xlabel('C_x (F)')
ylabel('P (W)')
legend('P_m_e_c_h, Rated Load 2 kW','P_m_e_c_h, Rated Load 4 kW','P_m_e_c_h, Rated Load 6 kW',...
    'P_l_o_a_d, Rated Load 2 kW','P_l_o_a_d, Rated Load 4 kW','P_l_o_a_d, Rated Load 6 kW')

figure(2)
plot(Cx_arr, f_cx(1,:),Cx_arr, f_cx(2,:),Cx_arr, f_cx(3,:))
title('External Excitation Capacitance vs. Frequency')
xlabel('C_x (F)')
ylabel('f (Hz)')
legend('Rated Load 2 kW','Rated Load 4 kW','Rated Load 6 kW')

figure(3)
plot3(Cx_arr,f_cx(1,:),Vph_cx(1,:),Cx_arr,f_cx(2,:),Vph_cx(2,:),Cx_arr,f_cx(3,:),Vph_cx(3,:))
title('Capacitance vs. Frequency vs. Phase Voltage')
grid
xlabel('C_x (F)')
ylabel('f (Hz)')
zlabel('V_p_h (V)')
legend('Rated Load 2 kW','Rated Load 4 kW','Rated Load 6 kW')

figure(4)
plot3(Cx_arr,-[Q_loss_cx(1,:).excitation], [Q_loss_cx(1,:).total], Cx_arr, -[Q_loss_cx(2,:).excitation], [Q_loss_cx(2,:).total], Cx_arr, -[Q_loss_cx(3,:).excitation], [Q_loss_cx(3,:).total])
title('Capacitance vs. Capacitive Q vs. Generator Q')
grid
xlabel('C_x (F)')
ylabel('Q_c_a_p (VAR)')
zlabel('Q_g_e_n (VAR)')
legend('Rated Load 2 kW','Rated Load 4 kW','Rated Load 6 kW')
ylim([0 1e4])
zlim([0 1e4])

%%
% % 3 phase electrical load
% P_load_in = 2e3:2e3:6e3;   %W
% Q_load = 0;         %VAR
% % rated values
% f_rated = 60;       %Hz
% Vline_rated = 220; %V
% poles = 2;
% % speeds
% % n_mech = 1.01*f_rated*120/poles;
n_mech_arr = linspace(2400,4000,500); %rpm
% % friction coefficients
% K_b = 0.4;
% K_w = 0.001;
% % stator values
% R1 = 0.148;         %ohm
% L1 = 0.423/2/pi/f_rated;     %H
% % rotor values
% R2 = 0.144;         %ohm
% L2 = 0.252/2/pi/f_rated;       %H
% % core values
% [magCurve.X, magCurve.Vg_f] = imageinterp;
% % Rm = 12000;         %ohm
% % Lm = 0.12;          %H
% % excitation values
% 
% Rx = 0;             %ohm
% Cx_arr = linspace(1e-6,400e-6,500);   %F
Cx = 240e-6;    %F

Vph_n = zeros(size(n_mech_arr));
f_n = zeros(size(n_mech_arr));
P_mech_n = zeros(size(n_mech_arr));
P_load_out_n = zeros(size(n_mech_arr));

% Test for different loads
for m = 1:length(P_load_in)
% Test varying capacitance
for n = 1:length(n_mech_arr);
% test with varying w but constant P_mech
[P_loss_n(m,n), Q_loss_n(m,n), Vph_n(m,n), f_n(m,n), P_mech_n(m,n), P_load_out_n(m,n) ] = SCIG_Ouazenne( P_load_in(m), Q_load, n_mech_arr(n), f_rated, Vline_rated, poles, K_b, K_w, R1, L1, R2, L2, magCurve, Rx, Cx );



end
end
unboundFlag = abs(Vph_n) > 1e3 | Vph_n < 0;
Vph_n(unboundFlag) = NaN;
f_n(unboundFlag) = NaN;
P_mech_n(unboundFlag) = NaN;
P_load_out_n(unboundFlag) = NaN;
% [Q_loss_n(unboundFlag).total] = [NaN];
% Q_loss_n(unboundFlag).excitation = NaN;

figure(5)
plot(   n_mech_arr, P_mech_n(1,:),n_mech_arr, P_mech_n(2,:),n_mech_arr, P_mech_n(3,:),...
        n_mech_arr, P_load_out_n(1,:),'-.',n_mech_arr, P_load_out_n(2,:),'-.',n_mech_arr, P_load_out_n(3,:),'-.')
title('Mech. Speed vs. Mech & Elec Power ')
xlabel('n (rpm)')
ylabel('P (W)')
legend('P_m_e_c_h, Rated Load 2 kW','P_m_e_c_h, Rated Load 4 kW','P_m_e_c_h, Rated Load 6 kW',...
    'P_l_o_a_d, Rated Load 2 kW','P_l_o_a_d, Rated Load 4 kW','P_l_o_a_d, Rated Load 6 kW')

figure(6)
plot(n_mech_arr, f_n(1,:),n_mech_arr, f_n(2,:),n_mech_arr, f_n(3,:))
title('Mechanical Speed vs. Frequency')
xlabel('n (rpm)')
ylabel('f (Hz)')
legend('Rated Load 2 kW','Rated Load 4 kW','Rated Load 6 kW')

figure(7)
plot3(n_mech_arr,f_n(1,:),Vph_n(1,:),n_mech_arr,f_n(2,:),Vph_n(2,:),n_mech_arr,f_n(3,:),Vph_n(3,:))
title('Mech. Speed vs. Frequency vs. Phase Voltage')
grid
xlabel('n (rpm)')
ylabel('f (Hz)')
zlabel('V_p_h (V)')
legend('Rated Load 2 kW','Rated Load 4 kW','Rated Load 6 kW')

figure(8)
plot(n_mech_arr,(f_n(1,:)/f_rated - n_mech_arr/3600)./(f_n(1,:)/f_rated), n_mech_arr, (f_n(2,:)/f_rated - n_mech_arr/3600)./(f_n(2,:)./f_rated), n_mech_arr, (f_n(3,:)/f_rated - n_mech_arr/3600)./(f_n(3,:)/f_rated))
title('Mech. Speed vs. Slip')
xlabel('n (rpm)')
ylabel('Slip')
legend('Rated Load 2 kW','Rated Load 4 kW','Rated Load 6 kW')
