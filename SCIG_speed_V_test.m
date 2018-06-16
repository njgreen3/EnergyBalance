clear 
close all
clc 

% This script tests the SCIG_Energy_Balance function under a range of 
% mechanical rotor speeds. The first test keeps the mechanical input power
% constant. The second test also varies the mechanical input power
% proportionally to the speed so torque is effectively kept constant.

P_mech = 1;
w_mech = .99:.001:1.15;

w_synch = 1;
K_b = .01;
K_w = .01;

R1 = .005;
X1 = .088;
R2 = .009;
X2 = .0125;
Rm = 1400;
Xm = 5;

% test with varying w but constant P 
[ P_out_p, Q_out_p, P_loss_p, Vph_p, E_p ] = SCIG_Energy_Balance( P_mech, w_mech, w_synch, K_b, K_w, R1, X1, R2, X2, Rm, Xm );
pf_p = P_out_p./sqrt(P_out_p.^2 + Q_out_p.^2);

figure(1)
plot(w_mech, abs(sqrt(3)*Vph_p),'b', w_mech, abs(sqrt(3)*E_p),'b--')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Voltage (per-unit)')

hold on

figure(2)
plot(w_mech, P_out_p,'b', w_mech, abs(Q_out_p),'b--')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Active and Reactive Power(per-unit)')

hold on

figure(3)
plot(w_mech,pf_p, 'b')
xlabel('Mechanical Anglular Velocity (per-unit)')
ylabel('Power Factor')

hold on

% test with varying w and proportionally varying P (constant torque)
P_mech = w_mech;
[ P_out_t, Q_out_t, P_loss_t, Vph_t, E_t ] = SCIG_Energy_Balance( P_mech, w_mech, w_synch, K_b, K_w, R1, X1, R2, X2, Rm, Xm );
pf_t = P_out_t./sqrt(P_out_t.^2 + Q_out_t.^2);

figure(1)
plot(w_mech, abs(sqrt(3)*Vph_t),'r', w_mech, sqrt(3)*abs(E_t),'r--')
legend('Phase voltage for constant mechanical input power', ...
    'Induced voltage for constant mechanical input power', ...
    'Phase voltage for constant input torque', ...
    'Induced voltage for constant input torque')

figure(2)
plot(w_mech, P_out_t,'r', w_mech, abs(Q_out_t),'r--')
legend('Active Power generated for constant mechanical input power', ...
    'Reactive Power consumed for constant mechanical input power', ...
    'Active Power generated for constant input torque', ...
    'Reactive Power consumed for constant input torque')

figure(3)
plot(w_mech, pf_t, 'r')
legend('Constant mechanical input power', ...
    'Constant input torque')