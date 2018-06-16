function [enthalpy] = pt2h(fluid_handle,pressure,temperature)
%PT2H Summary of this function goes here
%   This function takes caluculates the enthalpy of a fluid based off of
%   its temperature and pressure using the Cool Prob library. The fluid 
%   handle must already be generated and input in order to calculate the 
%   output values.

buffer_size = 255;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);
coder.extrinsic('calllib') 
coder.extrinsic('libpointer') 

[pt_ind,~] = calllib('coolprop','get_input_pair_index','PT_INPUTS'); %Mass specific enthalpy and pressure
[h_ind, ~] = calllib('coolprop','get_param_index','Hmass'); %temperature

var_len = length(pressure);
pressurePtr = libpointer('doublePtr',pressure);
temperaturePtr  = libpointer('doublePtr',temperature);
enthalpy = zeros(var_len,1);
enthaplyPtr = libpointer('doublePtr',enthalpy);

calllib('coolprop','AbstractState_update_and_1_out', fluid_handle, pt_ind, pressurePtr, temperaturePtr, var_len, h_ind, enthaplyPtr, ierr,herr,buffer_size);

enthalpy = get(enthaplyPtr,'Value');

end

