% VSI_test


P_rated = 1e4;
P_load = (0.1:0.01:1)*P_rated;

eff = InverterEfficiency(P_load, P_rated);
[P_source, P_loss] = Voltage_Source_Inverter_Simple(P_load,eff);

figure()
subplot(3,1,1)
plot(P_load, P_source)
ylabel('Source Power (W)')

subplot(3,1,2)
plot(P_load,P_loss)
ylabel('Power loss (W)')

subplot(3,1,3)
plot(P_load,eff)
ylabel('efficiency')

xlabel('Load Power (W)')