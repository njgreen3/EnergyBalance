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
fluid_wf = 'R245fa';
backend = 'HEOS';
[handle_source, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_source,ierr,herr,buffer_size);
[handle_sink, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_sink,ierr,herr,buffer_size);
[handle_wf, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_wf,ierr,herr,buffer_size);
