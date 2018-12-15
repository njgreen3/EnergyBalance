function [pt_plot_handle, hp_plot_handle, st_plot_handle] = R245faPTcurve(workingFluid)
% clear all

buffer_size = 1000;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);
backend = 'HEOS';

if nargin == 0
    workingFluid = 'R245fa';
%     workingFluid = 'R134a';
end

[handle, sh] = calllib('coolprop','AbstractState_factory',backend,workingFluid,ierr,herr,buffer_size);


% Temperature just above triple point is min temp for this curve
% T_min = ceil(CoolProp.PropsSI('Ttriple',' ',0,' ',0,workingFluid));
T_min = 270;

% Temperature just below critical point is max temp for this curve
% T_max = floor(CoolProp.PropsSI('Tcrit',' ',0,' ',0,workingFluid));
T_max = 380;

len = 1000;

% initialize temp and quality values and pointers for vaporization curve
T_in = linspace(T_min,T_max,len)';
Q_max = ones(size(T_in));
Q_min = zeros(size(T_in));

T_inPtr = libpointer('doublePtr',T_in);
Q_maxPtr = libpointer('doublePtr',Q_max);
Q_minPtr = libpointer('doublePtr',Q_min);
p_outPtr = libpointer('doublePtr',zeros(len,1));
H_maxPtr = libpointer('doublePtr',zeros(len,1));
H_minPtr = libpointer('doublePtr',zeros(len,1));
S_maxPtr = libpointer('doublePtr',zeros(len,1));
S_minPtr = libpointer('doublePtr',zeros(len,1));

% find input pair index for temp and quality
[qt_pair, ip] = calllib('coolprop','get_input_pair_index','QT_INPUTS');
% find parameter index for pressure
[p_ind, oi_p] = calllib('coolprop','get_param_index','P');
% find parameter index for enthalpy
[h_ind, oi_h] = calllib('coolprop','get_param_index','H');
% find parameter index for entropy
[s_ind, oi_h] = calllib('coolprop','get_param_index','S');


% initialize the temperature values for vaporization curve
% P_vap_curve = zeros(size(T_vap_curve));

% initialize the entropy values to 0
% S_curve = zeros(length(T_vap_curve));

% tic
calllib('coolprop','AbstractState_update_and_1_out', handle, qt_pair, Q_minPtr, T_inPtr, len, p_ind, p_outPtr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', handle, qt_pair, Q_maxPtr, T_inPtr, len, h_ind, H_maxPtr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', handle, qt_pair, Q_minPtr, T_inPtr, len, h_ind, H_minPtr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', handle, qt_pair, Q_maxPtr, T_inPtr, len, s_ind, S_maxPtr, ierr,herr,buffer_size);
calllib('coolprop','AbstractState_update_and_1_out', handle, qt_pair, Q_minPtr, T_inPtr, len, s_ind, S_minPtr, ierr,herr,buffer_size);
% toc
T_vap_curve = T_in;
P_vap_curve = get(p_outPtr,'Value');
H_vap_curve = get(H_maxPtr,'Value');
H_cond_curve = get(H_minPtr,'Value');
S_vap_curve = get(S_maxPtr,'Value');
S_cond_curve = get(S_minPtr,'Value');

% for n = 1:length(T_vap_curve)
%     P_vap_curve(n) = CoolProp.PropsSI('P','T',T_vap_curve(n),'Q',0,workingFluid);
% end
%
figure()
pt_plot_handle = subplot(3,1,1);
semilogy(T_vap_curve,P_vap_curve)
xlabel('Temperature (K)')
ylabel('Pressure (Pa)')
title(['p-T diagram for ' workingFluid ])

hp_plot_handle = subplot(3,1,2);
semilogy(H_vap_curve,P_vap_curve)
hold on
plot(H_cond_curve,P_vap_curve)
ylabel('Pressure (pa)')
xlabel('Mass Specific Enthalpy (J/kg)')
title(['p-h diagram for ' workingFluid ])
hold off

st_plot_handle = subplot(3,1,3);
plot(S_vap_curve,T_vap_curve)
hold on
plot(S_cond_curve,T_vap_curve)
ylabel('Temperature (K)')
xlabel('Mass Specific Entropy (J/kg-K)')
title(['S-T diagram for ' workingFluid ])
hold off

calllib('coolprop','AbstractState_free', handle, ierr,herr,buffer_size);
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