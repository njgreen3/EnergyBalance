% rotor speed
we = 2*pi*freq;
wr = we*linspace(.95,1.05,500); %range of speeds to generate Torque-Speed Plot
wrm = 2/pole_pairs*wr*(60/(2*pi)); %rotor speed (rpm)
s = (we - wr)/we; %Slip


% Thevenin voltage & impedance
Vth = 1i*Xm/(Rs + 1i*Xs + 1i*Xm)*Vline/sqrt(3);
Zth = 1/(1/(Rs + 1i*Xs) + 1/(1i*Xm));
Rth = real(Zth);
Xth = imag(Zth);

% torque and power vs. speed characteristics
Te = 3*(pole_pairs/2)*abs(Vth)^2*Rr./(s*we)./((Rth + Rr./s).^2 + (Xth + Xr)^2);
Pe = Te.*wr*(2/pole_pairs);

figure(1), clf
hold on
grid on
plot(wrm, Te)
xlabel('Rotor Speed, wrm (rpm)')
ylabel('Torque (Nm)')

figure(2)
plot(wrm, Pe, 'r')
hold on 
grid on
xlabel('Rotor Speed, wrm (rpm)')
ylabel('Electrical Power, (W)')