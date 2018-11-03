% This script plots pressure-temperature, pressure-enthalpy, and
% temperature-entropy 
[handle_temp, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_wf,ierr,herr,buffer_size);
[pt_handle, hp_handle, st_handle] = R245faPTcurve(fluid_wf);

% Values for verification
% case_ind = [24 49 74 100];
% cursor_shapecolor = ['xk'; 'ok'; 'sk'; 'dk';];

% Values for green field
% case_ind = [34 69 104 140];
% cursor_shapecolor = ['xk'; 'ok'; 'sk'; 'dk';];

% Values for brown field
% case_ind = [19 39 59 79];
% cursor_shapecolor = ['xk'; 'ok'; 'sk'; 'dk';];
case_ind = 50;
cursor_shapecolor = 'xk';

subplot(pt_handle)
hold on
for n = 1:length(case_ind)
     plot(T_wf.Data(case_ind(n),1:4),[p_hi.signals.values(2*n-1)*ones(1,2) p_low.signals.values(2*n-1)*ones(1,2)],cursor_shapecolor(n,:))
%     plot(T_wf.Data(case_ind(n),1:4),[p_hi.signals.values(1)*ones(1,2) p_low.signals.values(1)*ones(1,2)],cursor_shapecolor(n,:))

%     plot(T_wf.Data(case_ind(n),2),p_hi.signals.values(2*n),cursor_shapecolor(n,:))
%     plot(T_wf.Data(case_ind(n),3),p_low.signals.values(2*n),cursor_shapecolor(n,:))
%     plot(T_wf.Data(case_ind(n),4),p_low.signals.values(2*n),cursor_shapecolor(n,:))
%     plot(ORC_Temperature.pump_out,pressure_h,'xb')
end
hold off
% legend('Vaporization','Summer Low Flow','Summer High Flow','Winter Low Flow','Winter High Flow','Location','eastoutside')
legend('Vaporization','Brownfield Test')

subplot(hp_handle)
hold on
for n = 1:length(case_ind)

    plot(H_wf.Data(case_ind(n),1:4),[p_hi.signals.values(2*n-1)*ones(1,2) p_low.signals.values(2*n-1)*ones(1,2)],cursor_shapecolor(n,:))
%     plot(H_wf.Data(case_ind(n),1:4),[p_hi.signals.values(1)*ones(1,2) p_low.signals.values(1)*ones(1,2)],cursor_shapecolor(n,:))
%     plot(H_wf.Data(case_ind(n),2),p_hi.signals.values(2*n),cursor_shapecolor(n,:))
%     plot(H_wf.Data(case_ind(n),3),p_low.signals.values(2*n),cursor_shapecolor(n,:))
%     plot(H_wf.Data(case_ind(n),4),p_low.signals.values(2*n),cursor_shapecolor(n,:))
%     plot(H_pump_out,pressure_h,'xb')
end
hold off
% legend('Vaporization','Condensation','Summer Low Flow','Summer High Flow','Winter Low Flow','Winter High Flow','Location','eastoutside')
legend('Vaporization','Condensation','Brownfield Test')

% find parameter index for entropy
s_ind = calllib('coolprop','get_param_index','S');
% find input pair index for enthalpy and temperature
hp_pair = calllib('coolprop','get_input_pair_index','HmassP_INPUTS');

subplot(st_handle)
hold on
for n = 1:length(case_ind)
    H_temp = H_wf.Data(case_ind(n),:);
    H_ptr = libpointer('doublePtr',H_temp);
    p_temp = [p_hi.signals.values(2*n-1)*ones(1,2) p_low.signals.values(2*n-1)*ones(1,2)];
%     p_temp = [p_hi.signals.values(1)*ones(1,2) p_low.signals.values(1)*ones(1,2)];
    p_ptr = libpointer('doublePtr',p_temp);
    S_ptr = libpointer('doublePtr',zeros(size(H_temp)));
        

    
    calllib('coolprop','AbstractState_update_and_1_out', handle_temp, hp_pair, H_ptr, p_ptr, 4, s_ind, S_ptr, ierr,herr,buffer_size);
    S_temp = get(S_ptr,'Value');

    
    plot(S_temp(1:4),T_wf.Data(case_ind(n),1:4),cursor_shapecolor(n,:))
%     plot(S_temp(2),T_wf.Data(case_ind(n),2),cursor_shapecolor(n,:))
%     plot(S_temp(3),T_wf.Data(case_ind(n),3),cursor_shapecolor(n,:))
%     plot(S_temp(4),T_wf.Data(case_ind(n),4),cursor_shapecolor(n,:))
%     plot(S_pump_out,ORC_Temperature.pump_out,'xb')
end
hold off
% legend('Vaporization','Condensation','Summer Low Flow','Summer High Flow','Winter Low Flow','Winter High Flow','Location','eastoutside')
legend('Vaporization','Condensation','Brownfield Test')


calllib('coolprop','AbstractState_free', handle_temp, ierr,herr,buffer_size);
clear handle_temp T_temp H_temp T_ptr H_ptr