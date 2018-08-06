path_to_include= 'C:\Users\njgreen3\Documents\MATLAB\CoolPropStuff'; %specify path to coolprop's include folder

% Loading shared library
if ~libisloaded('coolprop') %checking whether library is already loaded
%     addpath(path_to_lib)
    addpath(path_to_include)
    libname = 'libCoolProp'; % OSX and linux
    if ispc
        libname = 'CoolProp';
    end
    loadlibrary(libname,'CoolPropLib.h','includepath',path_to_include,'alias','coolprop'); % loading library with alias coolprop
    disp('loaded CoolProp shared library.')
    disp('loaded these functions: ')
    libfunctions coolprop
end

% set error handling variables
buffer_size = 255;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);

% Close existing handles (to keep things tidy)
if exist('handle_source','var')
    calllib('coolprop','AbstractState_free', handle_source, ierr,herr,buffer_size);
end
if exist('handle_sink','var')
    calllib('coolprop','AbstractState_free', handle_sink, ierr,herr,buffer_size);
end
if exist('handle_wf','var')
    calllib('coolprop','AbstractState_free', handle_wf, ierr,herr,buffer_size);
end

% define the source, sink, and working fluids and set up respective handles
fluid_source = 'water';
fluid_sink = 'water';
fluid_wf = 'R134a';
backend = 'HEOS';
[handle_source, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_source,ierr,herr,buffer_size);
[handle_sink, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_sink,ierr,herr,buffer_size);
[handle_wf, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_wf,ierr,herr,buffer_size);

Tin_source = 346.48;
Tin_sink = 277.59;
m_dot_source = 33.39;
m_dot_sink = 101.68;
m_dot_wf_init = 12;
H_wf_init = 2.6e+05;
p_hi = 1695000;
p_low = 439000;
p_atm = 101325;

U_evap = 4500;
A_evap = 30;
U_cond = 5000;
A_cond = 50;

eff_turbine = 0.8;
eff_pump = 0.78;
eff_inverter = 0.93;
pf_grid = 0.9;

P_spd_speed = [834.591559026165	844.140855611678	853.690152197192	863.239448782706	872.788745368220	882.338041953733	891.887338539247	901.436635124761	910.985931710274	920.535228295788	930.084524881302	939.633821466816	949.183118052329	958.732414637843	968.281711223357	977.831007808871];
P_spd_power = -1/3*[166663.989372846	163937.628261044	155852.617467788	140628.948970298	116654.771019611	83123.6124212207	40806.0621676137	-7481.19463142690	-57239.4174183274	-103560.341265182	-142633.126036591	-172599.754334336	-193461.262910890	-206409.980979167	-213110.460435838	-215209.127126632];

Vline = 220;
freq = 60;
pole_pairs = 4;
Rs = 0.018;
Xs = 0.32;
Rr = 0.03;
Xr = 0.05;
Rm = 5000;
Xm = 18;
Rx = 0.003;
Xx = 12;


