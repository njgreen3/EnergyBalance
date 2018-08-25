% rotor speed
we = 2*pi*60;
wr = we'.*linspace(0,2,500); %range of speeds to generate Torque-Speed Plot
we = repmat(we',[1 500]);
wrm = 2/pole_pairs*wr*(60/(2*pi)); %rotor speed (rpm)
% wrm = 1/pole_pairs*wr*(60/(2*pi)); %rotor speed (rpm)
s = (we - wr)./we; %Slip


% Thevenin voltage & impedance
Vth = 1i*(Xm*2*pi./we)./(Rs + 1i*(Xs*2*pi./we) + 1i*(Xm*2*pi./we))*Vline/sqrt(3);
Zth = 1./(1./(Rs + 1i*(Xs*2*pi./we)) + 1./(1i*(Xm*2*pi./we)));
Rth = real(Zth);
Xth = imag(Zth);

% torque and power vs. speed characteristics
% Te = 3*(pole_pairs/1)*abs(Vth).^2*Rr./(s.*we)./((Rth + Rr./s).^2 + (Xth + Xr).^2);
% Pe = Te.*wr*(1/pole_pairs);
Te = 3*(pole_pairs/2)*abs(Vth).^2*Rr./(s.*we)./((Rth + Rr./s).^2 + (Xth + Xr).^2);
Pe = Te.*wr*(2/pole_pairs);


figure(1)
hold on
grid on
plot(wrm,Te)
xlabel('n (rpm)')
ylabel('Torque (Nm)')
% surfc(we/2/pi,wrm,Te,'EdgeColor','flat')
% xlabel('F_e (Hz)')
% ylabel('n (rpm)')
% zlabel('Torque (Nm)')

figure(2)
plot(wrm,Pe)
hold on 
grid on
xlabel('Rotor Speed, wrm (rpm)')
ylabel('Electrical Power, (W)')
% surfc(we/2/pi,wrm,Pe,'EdgeColor','flat')
% hold on 
% grid on
% xlabel('F_e (Hz)')
% ylabel('Rotor Speed, wrm (rpm)')
% zlabel('Electrical Power, (W)')