clear all
workingFluid = 'R245fa';

% Temperature just above triple point is min temp for this curve
% T_min = ceil(CoolProp.PropsSI('Ttriple',' ',0,' ',0,workingFluid));
T_min = 288;

% Temperature just below critical point is max temp for this curve
% T_max = floor(CoolProp.PropsSI('Tcrit',' ',0,' ',0,workingFluid));
T_max = 360;

% initialize the temperature values for vaporization curve
T_vap_curve = linspace(T_min,T_max,1000);

% initialize the temperature values for vaporization curve
P_vap_curve = zeros(size(T_vap_curve));

% initialize the entropy values to 0
% S_curve = zeros(length(T_vap_curve));
tic
for n = 1:length(T_vap_curve)
    P_vap_curve(n) = CoolProp.PropsSI('P','T',T_vap_curve(n),'Q',0,workingFluid);
end
toc
figure()
semilogy(T_vap_curve,P_vap_curve)
xlabel('Temperature (K)')
ylabel('Pressure (Pa)')
title(['ORC P-T diagram with ' workingFluid ' Vaporization Curve'])

% T_C = T_vap_curve - 273;
% P_atm = P_vap_curve * 9.86923e-6;
% 
% figure()
% semilogy(T_C,P_atm)
% xlabel('Temperature (C)')
% ylabel('Pressure (atm)')
% title(['ORC P-T diagram with ' workingFluid ' Vaporization Curve'])
% 
% 
% for n = 1:length(T_vap_curve)
%     for m = 1:length(T_vap_curve)
%         if n ~= m
%             S_curve(n,m) = CoolProp.PropsSI('S','T',T_vap_curve(n),'P',P_vap_curve(m),workingFluid);
%         end
%     end
% end
% 
% figure(2)
% surf(T_vap_curve,P_vap_curve,S_curve)

% pl = 1e5;
% ph = 1e6;
% 
% S_min = CoolProp.PropsSI('S','T',T_min,'P',pl,workingFluid);
% S_max = CoolProp.PropsSI('S','T',T_max,'P',ph,workingFluid);
% S_curve = linspace(S_min,S_max,length(T_vap_curve));
% T_pl_s_curve = zeros(length(T_vap_curve));
% T_ph_s_curve = zeros(length(T_vap_curve));
% 
% H_min = CoolProp.PropsSI('H','T',T_min,'P',pl,workingFluid);
% H_max = CoolProp.PropsSI('H','T',T_max,'P',ph,workingFluid);
% H_curve = linspace(H_min,H_max,length(T_vap_curve));
% T_pl_h_curve = zeros(length(T_vap_curve));
% T_ph_h_curve = zeros(length(T_vap_curve));
% 
% for n = 1:length(T_vap_curve)
%     T_pl_s_curve(n) = CoolProp.PropsSI('T','S',S_curve(n),'P',pl,workingFluid);
%     T_ph_s_curve(n) = CoolProp.PropsSI('T','S',S_curve(n),'P',ph,workingFluid);
%     
%     T_pl_h_curve(n) = CoolProp.PropsSI('T','H',H_curve(n),'P',pl,workingFluid);
%     T_ph_h_curve(n) = CoolProp.PropsSI('T','H',H_curve(n),'P',ph,workingFluid);
% end
%     