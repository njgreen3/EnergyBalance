function [temperature] = hp2t(fluid_handle,enthalpy,pressure)
%HP2T Summary of this function goes here
%   This function takes caluculates the temperature of a fluid based off of
%   its enthalpy and pressure using the Cool Prob library. The fluid handle
%   must already be generated and input in order to calculate the output
%   values.

buffer_size = 255;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);
coder.extrinsic('calllib') 
coder.extrinsic('libpointer') 

[hp_ind,~] = calllib('coolprop','get_input_pair_index','HmassP_INPUTS'); %Mass specific enthalpy and pressure
[t_ind, ~] = calllib('coolprop','get_param_index','T'); %temperature

var_len = length(pressure);
pressurePtr = libpointer('doublePtr',pressure);
enthaplyPtr = libpointer('doublePtr',enthalpy);
temperature = zeros(var_len,1);
temperaturePtr = libpointer('doublePtr',temperature);

calllib('coolprop','AbstractState_update_and_1_out', fluid_handle, hp_ind, enthaplyPtr, pressurePtr, var_len, t_ind, temperaturePtr, ierr,herr,buffer_size);

temperature = get(temperaturePtr,'Value');

end

