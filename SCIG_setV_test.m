clear 
close all
clc 
%%
% This script tests the SCIG_Energy_Balance function under a range of 
% mechanical rotor speeds. 

% electrical load
P_load = 0;
Q_load = 0;
Vph = 1;
% speeds
w_mech = 1.01;
w_synch = 1;
% friction coefficients
K_b = 0.01;
K_w = 0.01;
% stator values
R1 = 0.005;
L1 = 0.088;
% rotor values
R2 = 0.009;
L2 = 0.0125;
% core values
Rm = 1400;
Lm = 5;
% excitation values
Rx = 0;
Cx = 0;
% Cx = repmat(Cx,1,length(w_mech));


%%

w_mech = .99:.001:1.15;
% test with varying w but constant P_mech
[ P_loss_w, Q_loss_w, P_mech_w, S_grid_w ] = SCIG_Energy_Balance_setV( P_load, Q_load, Vph, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );


figure(1)
plot(w_mech, P_mech_w)
title('Mechanical Power vs. Rotor Velocity')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Power (per-unit)')



figure(2)
subplot(2,1,1)
plot(w_mech, P_loss_w.total)
title('Real and Reactive Power losses vs. Rotor Velocity')
ylabel('P_l_o_s_s (per-unit)')

subplot(2,1,2)
plot(w_mech, Q_loss_w.total)
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Q_l_o_s_s (per-unit)')

figure(3)
subplot(2,1,1)
plot(w_mech, real(S_grid_w))
title('Power Input vs. Rotor Velocity')
ylabel('Active Power (per-unit)')

subplot(2,1,2)
plot(w_mech, imag(S_grid_w))
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Reactive Power (per-unit)')



%%
% test with varying Cx
P_mech_w = 1;
w_mech = 1.01;
Cx = linspace(0,5,20);
% Cx = linspace(0,5,2000);
% Cx(Cx<.3 & Cx>.1) = [];

% test with varying Cx
[ P_loss_cx, Q_loss_cx, P_mech_cx, S_grid_cx ] = SCIG_Energy_Balance_setV( P_load, Q_load, Vph, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );

figure(4)
plot(Cx, P_mech_cx)
title('Mechanical Power vs. External Excitation Capacitance')
xlabel('Capacitance (per-unit)')
ylabel('Power (per-unit)')



figure(5)
subplot(2,1,1)
plot(Cx, P_loss_cx.total)
title('Real and Reactive Power losses vs. External Excitation Capacitance')
ylabel('P_l_o_s_s (per-unit)')

subplot(2,1,2)
plot(Cx, Q_loss_cx.total)
xlabel('Capacitance (per-unit)')
ylabel('Q_l_o_s_s (per-unit)')

figure(6)
subplot(2,1,1)
plot(Cx, real(S_grid_cx))
title('Power Input vs. External Excitation Capacitance')
ylabel('P_g_r_i_d (per-unit)')

subplot(2,1,2)
plot(Cx, imag(S_grid_cx))
xlabel('Capacitance (per-unit)')
ylabel('Q_g_r_i_d (per-unit)')


%%
% test with varying phase voltage
Cx = 0;
Vph = linspace(.8,1.2,501);


[ P_loss_v, Q_loss_v, P_mech_v, S_grid_v ] = SCIG_Energy_Balance_setV( P_load, Q_load, Vph, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );

figure(7)
plot(Vph, P_mech_v)
title('Mechanical Power vs. Phase Voltage')
xlabel('Voltage (per-unit)')
ylabel('Power (per-unit)')



figure(8)
subplot(2,1,1)
plot(Vph, P_loss_v.total)
title('Real and Reactive Power losses vs. Phase Voltage')
ylabel('P_l_o_s_s (per-unit)')

subplot(2,1,2)
plot(Vph, Q_loss_v.total)
xlabel('Voltage (per-unit)')
ylabel('Q_l_o_s_s (per-unit)')

figure(9)
subplot(2,1,1)
plot(Vph, real(S_grid_v))
title('Power Input vs. Phase Voltage')
ylabel('Active Power (per-unit)')

subplot(2,1,2)
plot(Vph, imag(S_grid_v))
xlabel('Voltage (per-unit)')
ylabel('Reactive Power (per-unit)')


%%
% 3d surface plots varying both speed and voltage
Vph = linspace(0.8,1.2,501);
w_mech = 1:.0002:1.035;
Vph3d = repmat(Vph', size(w_mech));
w_mech3d = repmat(w_mech, size(Vph'));

[ P_loss_v_w, Q_loss_v_w, P_mech_v_w, S_v_w] = SCIG_Energy_Balance_setV( P_load, Q_load, Vph3d, w_mech3d, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );

constP = 1;
constP_err = 0.002;
v_P1_0 = Vph3d(real(S_v_w)> (-constP - constP_err) & real(S_v_w) < (-constP + constP_err));
w_P1_0 = w_mech3d(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< (-constP + constP_err));
P_P1_0 = real(S_v_w(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< (-constP + constP_err)));
Q_P1_0 = imag(S_v_w(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< (-constP + constP_err)));
C_P1_0 = Q_P1_0./(w_P1_0.*v_P1_0.^2);

constP = 1.1;
v_P1_1 = Vph3d(real(S_v_w)> (-constP - constP_err) & real(S_v_w) < (-constP + constP_err));
w_P1_1 = w_mech3d(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< (-constP + constP_err));
P_P1_1 = real(S_v_w(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< (-constP + constP_err)));
Q_P1_1 = imag(S_v_w(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< (-constP + constP_err)));
C_P1_1 = Q_P1_1./(w_P1_1.*v_P1_1.^2);

constP = 0.9;
v_P0_9 = Vph3d(real(S_v_w)> (-constP - constP_err) & real(S_v_w) < -constP + constP_err);
w_P0_9 = w_mech3d(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< -constP + constP_err);
P_P0_9 = real(S_v_w(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< -constP + constP_err));
Q_P0_9 = imag(S_v_w(real(S_v_w)> (-constP - constP_err) & real(S_v_w)< -constP + constP_err));
C_P0_9 = Q_P0_9./(w_P0_9.*v_P0_9.^2);

figure(10)
surfc(Vph3d, w_mech3d, P_mech_v_w,'EdgeColor','flat')
title('Mechanical Power vs Phase Voltage vs. Rotational Speed')
zlabel('Power (Per-Unit)')
ylabel('Mechanical Anglular Velocity (per-unit)')
xlabel('Voltage (per-unit)')

figure(11)
subplot(2,1,1)
surfc(Vph3d, w_mech3d, P_loss_v_w.total,'EdgeColor','flat')
title('Power Loss vs Phase Voltage vs. Rotational Speed')
zlabel('Active Power (Per-Unit)')
ylabel('Mechanical Anglular Velocity (per-unit)')
xlabel('Voltage (per-unit)')
subplot(2,1,2)
surfc(Vph3d, w_mech3d, Q_loss_v_w.total,'EdgeColor','flat')
zlabel('Reactive Power (Per-Unit)')
ylabel('Mechanical Anglular Velocity (per-unit)')
xlabel('Voltage (per-unit)')

figure(12)
subplot(2,1,1)
surfc(Vph3d, w_mech3d, real(S_v_w),'EdgeColor','flat')
title('Input Power vs Phase Voltage vs. Rotational Speed')
zlabel('Active Power (Per-Unit)')
ylabel('Mechanical Anglular Velocity (per-unit)')
xlabel('Voltage (per-unit)')
subplot(2,1,2)
surfc(Vph3d, w_mech3d, imag(S_v_w),'EdgeColor','flat')
zlabel('Reactive Power (Per-Unit)')
ylabel('Mechanical Anglular Velocity (per-unit)')
xlabel('Voltage (per-unit)')

figure(13)
plot(v_P0_9,w_P0_9, v_P1_0,w_P1_0, v_P1_1,w_P1_1)
title('Speed vs Voltage under constant output power')
xlabel('Voltage (per-unit)')
ylabel('Mechanical Anglular Velocity (per-unit)')
legend('P = 0.9 pu','P = 1.0 pu','P = 1.1 pu')

figure(14)
plot(v_P0_9,C_P0_9, v_P1_0,C_P1_0, v_P1_1,C_P1_1)
title('External Capacitance vs Voltage under constant output power')
xlabel('Voltage (per-unit)')
ylabel('Excitation Capacitance (per-unit)')
legend('P = 0.9 pu','P = 1.0 pu','P = 1.1 pu')

figure(15)
plot(v_P0_9,Q_P0_9, v_P1_0,Q_P1_0, v_P1_1,Q_P1_1)
title('External Reactive Power vs Voltage under constant output power')
xlabel('Voltage (per-unit)')
ylabel('Reactive Power (per-unit)')
legend('P = 0.9 pu','P = 1.0 pu','P = 1.1 pu')

%%
% 3d surface plots varying both Cx and voltage
w_mech = 1.01;
Vph = linspace(.8,1.2,501);
Cx = linspace(0,1,501);
Vph3d = repmat(Vph', size(Cx));
Cx3d = repmat(Cx, size(Vph'));

[ P_loss_v_c, Q_loss_v_c, P_mech_v_c, S_v_c] = SCIG_Energy_Balance_setV( P_load, Q_load, Vph3d, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx3d );

figure(16)
surfc(Vph3d, Cx3d, P_mech_v_c,'EdgeColor','flat')
title('Mechanical Power vs Phase Voltage vs. External Capacitance')
zlabel('Power (Per-Unit)')
ylabel('Capacitance (per-unit)')
xlabel('Voltage (per-unit)')

figure(17)
subplot(2,1,1)
surfc(Vph3d, Cx3d, P_loss_v_c.total,'EdgeColor','flat')
title('Power Loss vs Phase Voltage vs. External Capacitance')
zlabel('Active Power (Per-Unit)')
ylabel('Capacitance (per-unit)')
xlabel('Voltage (per-unit)')
subplot(2,1,2)
surfc(Vph3d, Cx3d, Q_loss_v_c.total,'EdgeColor','flat')
zlabel('Reactive Power (Per-Unit)')
ylabel('Capacitance (per-unit)')
xlabel('Voltage (per-unit)')

figure(18)
subplot(2,1,1)
surfc(Vph3d, Cx3d, real(S_v_c),'EdgeColor','flat')
title('Input Reactive Power vs Phase Voltage vs. External Capacitance')
zlabel('Active Power (Per-Unit)')
ylabel('Capacitance (per-unit)')
xlabel('Voltage (per-unit)')
subplot(2,1,2)
surfc(Vph3d, Cx3d, imag(S_v_c),'EdgeColor','flat')
zlabel('Reactive Power (Per-Unit)')
ylabel('Capacitance (per-unit)')
xlabel('Voltage (per-unit)')
