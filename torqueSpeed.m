% rotor speed
we = 2*pi*(30:5:70);
wr = we'.*linspace(.9,1.1,500); %range of speeds to generate Torque-Speed Plot
we = repmat(we',[1 500]);
wrm = 2/pole_pairs*wr*(60/(2*pi)); %rotor speed (rpm)
s = (we - wr)./we; %Slip


% Thevenin voltage & impedance
Vth = 1i*(Xm*2*pi./we)./(Rs + 1i*(Xs*2*pi./we) + 1i*(Xm*2*pi./we))*Vline/sqrt(3);
Zth = 1./(1./(Rs + 1i*(Xs*2*pi./we)) + 1./(1i*(Xm*2*pi./we)));
Rth = real(Zth);
Xth = imag(Zth);

% torque and power vs. speed characteristics
Te = 3*(pole_pairs/2)*abs(Vth).^2*Rr./(s.*we)./((Rth + Rr./s).^2 + (Xth + Xr).^2);
Pe = Te.*wr*(2/pole_pairs);

figure(1), clf
hold on
grid on
surfc(we/2/pi,wrm,Te,'EdgeColor','flat')
xlabel('F_e (Hz)')
ylabel('n (rpm)')
zlabel('Torque (Nm)')

figure(2)
surfc(we/2/pi,wrm,Pe,'EdgeColor','flat')
hold on 
grid on
xlabel('F_e (Hz)')
ylabel('Rotor Speed, wrm (rpm)')
zlabel('Electrical Power, (W)')