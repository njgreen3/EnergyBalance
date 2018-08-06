% clear all

buffer_size = 1000;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);
backend = 'HEOS';

% workingFluid = 'R245fa';
workingFluid = 'R134a';

[handle, sh] = calllib('coolprop','AbstractState_factory',backend,workingFluid,ierr,herr,buffer_size);


% Temperature just above triple point is min temp for this curve
% T_min = ceil(CoolProp.PropsSI('Ttriple',' ',0,' ',0,workingFluid));
T_min = 280;

% Temperature just below critical point is max temp for this curve
% T_max = floor(CoolProp.PropsSI('Tcrit',' ',0,' ',0,workingFluid));
T_max = 360;

len = 1000;

% initialize temp and quality values and pointers for vaporization curve
T_in = linspace(T_min,T_max,len)';
Q_in = ones(size(T_in));

T_inPtr = libpointer('doublePtr',T_in);
Q_inPtr = libpointer('doublePtr',Q_in);
p_outPtr = libpointer('doublePtr',zeros(len,1));

% find input pair index for temp and quality
[input_pair, ip] = calllib('coolprop','get_input_pair_index','QT_INPUTS');
% find parameter index for pressure
[output_ind, oi] = calllib('coolprop','get_param_index','P');



% initialize the temperature values for vaporization curve
% P_vap_curve = zeros(size(T_vap_curve));

% initialize the entropy values to 0
% S_curve = zeros(length(T_vap_curve));

tic
calllib('coolprop','AbstractState_update_and_1_out', handle, input_pair, Q_inPtr, T_inPtr, len, output_ind, p_outPtr, ierr,herr,buffer_size);
toc
T_vap_curve = T_in;
P_vap_curve = get(p_outPtr,'Value');

% for n = 1:length(T_vap_curve)
%     P_vap_curve(n) = CoolProp.PropsSI('P','T',T_vap_curve(n),'Q',0,workingFluid);
% end
% 
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