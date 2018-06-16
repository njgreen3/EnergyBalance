clear 
close all
clc 
%%
% This script tests the SCIG_Admittance_Balance function under a range of 
% inputs inculding external excitation capacitance, rotor speeds, and
% loads.

% 3 phase electrical load
P_load = [4e6;4.5e6;5e6];   %W
Q_load = 0;         %VAR
% rated values
f_rated = 60;       %Hz
Vline_rated = 6600; %V
poles = 6;
% speeds
n_mech = 1.01*f_rated*120/poles;
% friction coefficients
K_b = 0.5;
K_w = 0.05;
% stator values
R1 = 0.044;         %ohm
L1 = 0.002;         %H
% rotor values
R2 = 0.078;         %ohm
L2 = 0.00029;       %H
% core values
Rm = 12000;         %ohm
Lm = 0.12;          %H
% excitation values
Rx = 0;             %ohm
Cx_arr = linspace(5e-5,5e-4,200);   %F
% Cx = repmat(Cx,1,length(w_mech));

Vph = zeros(size(Cx_arr));
f = zeros(size(Cx_arr));
P_mech = zeros(size(Cx_arr));

% Test for different loads
for m = 1:length(P_load)
% Test varying capacitance
for n = 1:length(Cx_arr);
% test with varying w but constant P_mech
[P_loss(n,m), Q_loss(n,m), Vph_cx(n,m), f_cx(n,m), P_mech_cx(n,m) ] = SCIG_Admittance_Balance( P_load(m), Q_load, n_mech, f_rated, Vline_rated, poles, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx_arr(n) );

end
end

figure(1)
plot(Cx_arr, P_mech_cx(:,1),Cx_arr, P_mech_cx(:,2),Cx_arr, P_mech_cx(:,3))
title('External Excitation Capacitance vs. Mechanical Power')
xlabel('C_x (F)')
ylabel('P_m_e_c_h (V)')
legend('P = 4 MW','P = 4.5 MW','P = 5 MW')

figure(2)
plot(Cx_arr, f_cx(:,1),Cx_arr, f_cx(:,2),Cx_arr, f_cx(:,3))
title('External Excitation Capacitance vs. Frequency')
xlabel('C_x (F)')
ylabel('f (Hz)')
legend('P = 4 MW','P = 4.5 MW','P = 5 MW')

figure(3)
plot3(Cx_arr,f_cx(:,1),Vph_cx(:,1),Cx_arr,f_cx(:,2),Vph_cx(:,2),Cx_arr,f_cx(:,3),Vph_cx(:,3))
title('Capacitance vs. Frequency vs. Phase Voltage')
grid
xlabel('C_x (F)')
ylabel('f (Hz)')
zlabel('V_p_h (V)')
legend('P = 4 MW','P = 4.5 MW','P = 5 MW')

figure(4)
plot3(Cx_arr,-[Q_loss(:,1).excitation], [Q_loss(:,1).total],Cx_arr,-[Q_loss(:,2).excitation], [Q_loss(:,2).total],Cx_arr,-[Q_loss(:,3).excitation], [Q_loss(:,3).total])
title('Capacitance vs. Capacitive Q vs. Generator Q')
grid
xlabel('C_x (F)')
ylabel('Q_c (F)')
zlabel('Q_g_e_n (F)')
legend('P = 4 MW','P = 4.5 MW','P = 5 MW')
% %%
% 
% % 3 phase electrical load
% P_load = 4e6;       %W
% Q_load = 0;         %VAR
% % rated values
% f_rated = 60;       %Hz
% Vline_rated = 6600; %V
% poles = 6;
% % speeds
% n_mech_arr = f_rated*120/poles*linspace(1,1.09,200);    %rpm
% % friction coefficients
% K_b = 0.5;
% K_w = 0.05;
% % stator values
% R1 = 0.044;         %ohm
% L1 = 0.002;         %H
% % rotor values
% R2 = 0.078;         %ohm
% L2 = 0.00029;       %H
% % core values
% Rm = 12000;         %ohm
% Lm = 0.12;          %H
% % excitation values
% Rx = 0;             %ohm
% Cx = 7.7e-5;        %F
% % Cx = repmat(Cx,1,length(w_mech));
% 
% Vph_nmech = zeros(size(n_mech_arr));
% f_nmech = zeros(size(n_mech_arr));
% P_mech_nmech = zeros(size(n_mech_arr));
% 
% % Test varying mech speed
% for n = 1:length(n_mech_arr);
% % test with varying w but constant P_mech
% [P_loss(n), Q_loss(n), Vph_nmech(n), f_nmech(n), P_mech_nmech(n) ] = SCIG_Admittance_Balance( P_load, Q_load, n_mech_arr(n), f_rated, Vline_rated, poles, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );
% 
% end
% 
% figure(5)
% plot(n_mech_arr, P_mech_nmech)
% title('Mechanical Speed vs. Mechanical Power')
% xlabel('n_m_e_c_h (rpm)')
% ylabel('P_m_e_c_h (V)')
% 
% figure(6)
% plot(n_mech_arr, f_nmech)
% title('Mechanical Speed vs. Frequency')
% xlabel('n_m_e_c_h (rpm)')
% ylabel('f (Hz)')
% 
% figure(7)
% plot3(n_mech_arr,f_nmech,Vph_nmech)
% title('Mechanical Speed vs. Frequency vs. Phase Voltage')
% grid
% xlabel('n_m_e_c_h (rpm)')
% ylabel('f (Hz)')
% zlabel('V_p_h (V)')
% 
% figure(8)
% plot3(n_mech_arr,-[Q_loss.excitation], [Q_loss.total])
% title('Mechanical Speed vs. Capacitive Q vs. Generator Q')
% grid
% xlabel('n_m_e_c_h (rpm)')
% ylabel('Q_c (F)')
% zlabel('Q_g_e_n (F)')