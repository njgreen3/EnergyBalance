clear 
close all
clc 
%%
% This script tests the SCIG_Energy_Balance function under a range of 
% mechanical rotor speeds. The first test keeps the mechanical input power
% constant. The second test also varies the mechanical input power
% proportionally to the speed so torque is effectively kept constant.

P_mech = 1;

% electrical load
P_load = 0.9;
Q_load = 0;
% speeds
% w_mech = .99:.01:1.15;
w_mech = .99:.001:1.15;
% w_mech(w_mech == 1) = [];
% w_mech(2:5,:) = [w_mech;w_mech;w_mech;w_mech];
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
Cx = 1/.74;
% Cx = repmat(Cx,1,length(w_mech));


% % test with varying w but constant P_mech
% [P_loss(n), Q_loss(n), Vph(n), f(n), P_mech(n) ] = SCIG_Admittance_Balance( P_load, Q_load, n_mech, f_rated, Vline_rated, poles, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx(n) );
% 
% 
% figure(1)
% plot(w_mech, abs(Vph_p),'b', w_mech, abs(E_p),'b--')
% title('Voltage vs. Rotor Velocity')
% xlabel('Mechanical Anglular Velocity (per-unit)')
% ylabel('Voltage (per-unit)')
% 
% hold on
% 
% figure(2)
% subplot(2,1,1)
% plot(w_mech, P_loss_p.total)
% title('Real and Reactive Power losses vs. Rotor Velocity')
% ylabel('P_l_o_s_s (per-unit)')
% hold on
% 
% subplot(2,1,2)
% plot(w_mech, Q_loss_p.total)
% xlabel('Mechanical Anglular Velocity (per-unit)')
% ylabel('Q_l_o_s_s (per-unit)')
% hold on

% test with varying w and proportionally varying P_mech (constant torque)
P_mech = w_mech;
[ P_loss_t, Q_loss_t, Vph_t, E_t ] = SCIG_Energy_Balance_excitation( P_mech, P_load, Q_load, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );


figure(1)
plot(w_mech, abs(Vph_t),'r', w_mech, abs(E_t),'r--')
legend('|V_p_h| constant mechanical input power', ...
    '|E| constant mechanical input power', ...
    '|V_p_h| constant input torque', ...
    '|E| constant input torque')

figure(2)
subplot(2,1,1)
plot(w_mech, P_loss_t.total)
legend('Constant mechanical input power', ...
    'Constant input torque')

subplot(2,1,2)
plot(w_mech, Q_loss_t.total)
legend('Constant mechanical input power', ...
    'Constant input torque')


%%
% test with varying Cx
P_mech = 1;
w_mech = 1.04;
% Cx = linspace(1,5,200);
Cx = linspace(0,5,2000);
% Cx(Cx<.3 & Cx>.1) = [];

% test with varying Cx
[ P_loss_cx, Q_loss_cx, Vph_cx, E_cx ] = SCIG_Energy_Balance_excitation( P_mech, P_load, Q_load, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );

figure(3)
plot(Cx, abs(Vph_cx), Cx, abs(E_cx))
title('Voltage vs. Excitation Shunt Capacitance')
xlabel('Capacitance (per-unit)')
ylabel('Voltage (per-unit)')
legend('|V_p_h|', '|E|')
% ylim([0 2])
% xlim([0 5])
hold on

figure(4)
subplot(2,1,1)
plot(Cx, P_loss_cx.total)
title('Real and Reactive Power losses vs. Excitation Shunt Capacitance')
ylabel('P_l_o_s_s (per-unit)')
ylim([0 1])
hold on
subplot(2,1,2)
plot(Cx, Q_loss_cx.total)
ylabel('Q_l_o_s_s (per-unit)')
xlabel('Capacitance (per-unit)')
% legend('P_l_o_s_s', 'Q_l_o_s_s')
ylim([0 2])
hold on

figure(5)
plot(1./(Cx*w_synch), abs(Vph_cx), 1./(Cx.*w_synch), abs(E_cx))
title('Voltage vs. Excitation Shunt Impedance')
xlabel('Reactive Impedance (per-unit)')
ylabel('Voltage (per-unit)')
legend('|V_p_h|', '|E|')
% xlim([0 2])
hold on

figure(6)
subplot(2,1,1)
plot(1./(Cx.*w_synch), P_loss_cx.total)
title('Real and Reactive Power losses vs. Excitation Shunt Impedance')
ylabel('P_l_o_s_s (per-unit)')
% xlim([0 2])
hold on
subplot(2,1,2)
plot(1./(Cx.*w_synch), Q_loss_cx.total)
ylabel('Q_l_o_s_s (per-unit)')
xlabel('Reactive Impedance (per-unit)')
% legend('P_l_o_s_s ', 'Q_l_o_s_s ')
% xlim([0 2])
% ylim([-4 2])
hold on


%%
% test with varying load impedance angle with constant S
Cx = 1/.74;

impedanceAngle = linspace(0,90,200);
S_load = 0.9;
% S_load = P_load;
P_load = S_load*cos(impedanceAngle*pi/180);
Q_load = S_load*sin(impedanceAngle*pi/180);

% test with varying load power angle
[ P_loss_pa, Q_loss_pa, Vph_pa, E_pa ] = SCIG_Energy_Balance_excitation( P_mech, P_load, Q_load, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );

figure(7)
plot(impedanceAngle, abs(Vph_pa), impedanceAngle, abs(E_pa))
title('Voltage vs. Load Impedance Angle')
xlabel('Impedance Angle (deg)')
ylabel('Voltage (per-unit)')
legend('|V_p_h|', '|E|')
% ylim([0 2])
hold on

figure(8)
subplot(2,1,1)
plot(impedanceAngle, P_loss_pa.total)
title('Real and Reactive Power losses vs. Load Impedance Angle')
ylabel('P_l_o_s_s (per-unit)')
hold on
subplot(2,1,2)
plot(impedanceAngle, Q_loss_pa.total)
xlabel('Impedance Angle (deg)')
ylabel('Q_l_o_s_s (per-unit)')
% legend('P_l_o_s_s ', 'Q_l_o_s_s ')
% ylim([0 5])
hold on

%%
% varying load power with constant power factor

pf = 1;
load_theta = acos(pf);
P_load = linspace(0,0.95,200);
Q_load = P_load*tan(load_theta);
S_load = sqrt(P_load.^2 + Q_load.^2);

[ P_loss_pl, Q_loss_pl, Vph_pl, E_pl ] = SCIG_Energy_Balance_excitation( P_mech, P_load, Q_load, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );

figure(9)
plot(P_load, abs(Vph_pl), P_load, abs(E_pl))
title('Voltage vs. Load Power ')
xlabel('Load Power (Per-Unit)')
ylabel('Voltage (per-unit)')
legend('|V_p_h|', '|E|')
% ylim([0 2])
hold on

figure(10)
subplot(2,1,1)
plot(P_load, P_loss_pl.total)
title('Real and Reactive Power losses vs. Load Power')
ylabel('P_l_o_s_s (per-unit)')
hold on
subplot(2,1,2)
plot(P_load, Q_loss_pl.total)
xlabel('Load Power (Per-Unit)')
ylabel('Q_l_o_s_s (per-unit)')
hold on

%%
% 3d surface plots varying both load power (with constant pf) and Cx
P_load = linspace(0,0.95,50);
Cx = linspace(1,5,200);
Cx3d = repmat(Cx', size(P_load));
P_load3d = repmat(P_load, size(Cx'));
% Q_load3d = P_load3d*tan(load_theta);
Q_load3d = zeros(size(P_load3d));
vMax = 2;
pMax = .25;
qMax = 10;

[ P_loss_cxpl, Q_loss_cxpl, Vph_cxpl, E_cxpl ] = SCIG_Energy_Balance_excitation( P_mech, P_load3d, Q_load3d, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx3d );
% Vph_cxpl(Vph_cxpl > vMax) = vMax;
% E_cxpl(E_cxpl > vMax) = vMax;
% P_loss_cxpl.total(P_loss_cxpl.total > pMax) = pMax;
% Q_loss_cxpl.total(Q_loss_cxpl.total > qMax) = qMax;

figure(11)
surfc(P_load3d, Cx3d, abs(Vph_cxpl),'EdgeColor','flat')
title('Phase Voltage vs. Excitation Shunt Capacitance vs. Active Load Power')
xlabel('Power (Per-Unit)')
ylabel('Capacitance')
zlabel('Voltage (per-unit)')
% zlim([0 vMax])
% legend('|V_p_h|', '|E|')

figure(12)
surf(P_load3d, Cx3d, abs(E_cxpl),'EdgeColor','flat')
title('Internal Voltage vs. Excitation Shunt Capacitance vs. Active Load Power')
xlabel('Power (Per-Unit)')
ylabel('Capacitance')
zlabel('Voltage (per-unit)')
% zlim([0 vMax])
% legend('|V_p_h|', '|E|')
% ylim([0 2])

figure(13)
subplot(2,1,1)
surf(P_load3d, Cx3d, P_loss_cxpl.total,'EdgeColor','flat')
title('Real and Reactive Power losses vs. Excitation Shunt Capacitance vs. Active Load Power')
xlabel('Load Power (Per-Unit)')
ylabel('Capacitance')
zlabel('P_l_o_s_s (per-unit)')
subplot(2,1,2)
surf(P_load3d, Cx3d, Q_loss_cxpl.total,'EdgeColor','flat')
xlabel('Load Power (Per-Unit)')
ylabel('Capacitance')
zlabel('Q_l_o_s_s (per-unit)')

%%
% 3d surface plots varying both load power (with constant pf) and load impedance angle
P_load = linspace(0,0.95,96);
impedanceAngle = linspace(0,50,42);
Cx = 1/.74;

impedanceAngle3d = repmat(impedanceAngle', size(P_load));
P_load3d = repmat(P_load, size(impedanceAngle'));
Q_load3d = P_load3d.*tan(impedanceAngle3d*pi/180);
S_load3d = sqrt(P_load3d.^2 + Q_load3d.^2);

vMax = 2;
pMax = .25;
qMax = 10;

[ P_loss_pla, Q_loss_pla, Vph_pla, E_pla ] = SCIG_Energy_Balance_excitation( P_mech, P_load3d, Q_load3d, w_mech, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx );
% Vph_cxpl(Vph_cxpl > vMax) = vMax;
% E_cxpl(E_cxpl > vMax) = vMax;
% P_loss_cxpl.total(P_loss_cxpl.total > pMax) = pMax;
% Q_loss_cxpl.total(Q_loss_cxpl.total > qMax) = qMax;

figure(14)
surf(P_load3d, impedanceAngle3d, abs(Vph_pla),'EdgeColor','flat')
title('Phase Voltage vs. Load Impedance Angle vs. Active Load Power')
xlabel('Power (Per-Unit)')
ylabel('Impedance Angle (deg)')
zlabel('Voltage (per-unit)')
% zlim([0 vMax])
% legend('|V_p_h|', '|E|')

figure(15)
surf(P_load3d, impedanceAngle3d, abs(E_pla),'EdgeColor','flat')
title('Internal Voltage vs. Load Impedance Angle vs. Active Load Power')
xlabel('Power (Per-Unit)')
ylabel('Impedance Angle (deg)')
zlabel('Voltage (per-unit)')
% zlim([0 vMax])
% legend('|V_p_h|', '|E|')
% ylim([0 2])

figure(16)
subplot(2,1,1)
surf(P_load3d, impedanceAngle3d, P_loss_pla.total,'EdgeColor','flat')
title('Real and Reactive Power losses vs. Load Impedance Angle vs. Active Load Power')
xlabel('Load Power (Per-Unit)')
ylabel('Impedance Angle (deg)')
zlabel('P_l_o_s_s (per-unit)')
subplot(2,1,2)
surf(P_load3d, impedanceAngle3d, Q_loss_pla.total,'EdgeColor','flat')
xlabel('Load Power (Per-Unit)')
ylabel('Impedance Angle (deg)')
zlabel('Q_l_o_s_s (per-unit)')

%%
% 3d surface plots varying both w_mech (with constant input power) and Cx
w_mech = .99:.001:1.15;
% Cx = linspace(1,5,200);
Cx = linspace(0,.1,200);
Cx3d = repmat(Cx', size(w_mech));
w_mech3d = repmat(w_mech, size(Cx'));

pf = 1;
load_theta = acos(pf);
P_load = 0.90;
Q_load = P_load*tan(load_theta);
vMax = 2;
pMax = .25;
qMax = 10;

[ P_loss_cxw, Q_loss_cxw, Vph_cxw, E_cxw ] = SCIG_Energy_Balance_excitation( P_mech, P_load, Q_load, w_mech3d, w_synch, K_b, K_w, R1, L1, R2, L2, Rm, Lm, Rx, Cx3d );
% Vph_cxpl(Vph_cxpl > vMax) = vMax;
% E_cxpl(E_cxpl > vMax) = vMax;
% P_loss_cxpl.total(P_loss_cxpl.total > pMax) = pMax;
% Q_loss_cxpl.total(Q_loss_cxpl.total > qMax) = qMax;

figure(17)
surf(w_mech3d, Cx3d, abs(Vph_cxw),'EdgeColor','flat')
title('Phase Voltage vs. Excitation Shunt Capacitance vs. Rotor Speed')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Capacitance')
zlabel('Voltage (per-unit)')
% zlim([0 vMax])
% legend('|V_p_h|', '|E|')

figure(18)
surf(w_mech3d, Cx3d, abs(E_cxw),'EdgeColor','flat')
title('Internal Voltage vs. Excitation Shunt Capacitance vs. Rotor Speed')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Capacitance')
zlabel('Voltage (per-unit)')
% zlim([0 vMax])
% legend('|V_p_h|', '|E|')
% ylim([0 2])

figure(19)
subplot(2,1,1)
surf(w_mech3d, Cx3d, P_loss_cxw.total,'EdgeColor','flat')
title('Real and Reactive Power losses vs. Excitation Shunt Capacitance vs. Rotor Speed')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Capacitance')
zlabel('P_l_o_s_s (per-unit)')
subplot(2,1,2)
surf(w_mech3d, Cx3d, Q_loss_cxw.total,'EdgeColor','flat')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Capacitance')
zlabel('Q_l_o_s_s (per-unit)')