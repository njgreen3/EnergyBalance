function [ unit_val ] = pu2unit( desired_unit, pu_val, base_unit, base_val, varargin)
%PU2UNIT Summary of this function goes here
%   Detailed explanation goes here

% If there is an odd number of input arguments then something is wrong
if mod(nargin,2)
    error('Each input value must also have corresponding name.')
end

powerNames = {'kw','w','kvar','var','kva','va','power','p','q','s','watt','volt amp','volt amp reactive'};
voltageNames = {'v','kv','volt'};
currentNames = {'i','a','current','amp'};
impedanceNames = {'ohm','r','x','z','impedance','reactance','resistance'};
frequencyNames = {'hz','hertz','frequency','freq','f'};
capacitanceNames = {'farad','c','capacitance'};
inductanceNames = {'l','h','henry','inductance'};

baseFlag = 0;

varargin{end + 1} = base_unit;
varargin{end + 1} = base_val;

% check the name for each value in varagin
for n = 1:2:length(varargin)
    switch lower(varargin{n})
        case powerNames
            if bitget(baseFlag,7)
                warning('Multiple Power bases in input. Using first instance.')
                continue
            end
            s_base = varargin{n+1};
            baseFlag = bitset(baseFlag,7);
        case voltageNames
            if bitget(baseFlag,6)
                warning('Multiple Voltage bases in input. Using first instance.')
                continue
            end
            v_base = varargin{n+1};
            baseFlag = bitset(baseFlag,6);
        case currentNames
            if bitget(baseFlag,5)
                warning('Multiple Current bases in input. Using first instance.')
                continue
            end
            i_base = varargin{n+1};
            baseFlag = bitset(baseFlag,5);
        case impedanceNames
            if bitget(baseFlag,4)
                warning('Multiple Impedance bases in input. Using first instance.')
                continue
            end
            z_base = varargin{n+1};
            baseFlag = bitset(baseFlag,4);
        case frequencyNames
            if bitget(baseFlag,3)
                warning('Multiple Frequency bases in input. Using first instance.')
                continue
            end
            f_base = varargin{n+1};
            baseFlag = bitset(baseFlag,3);
        case capacitanceNames
            if bitget(baseFlag,2)
                warning('Multiple Capacitance bases in input. Using first instance.')
                continue
            end
            c_base = varargin{n+1};
            baseFlag = bitset(baseFlag,2);
        case inductanceNames
            if bitget(baseFlag,1)
                warning('Multiple Inductance bases in input. Using first instance.')
                continue
            end
            l_base = varargin{n+1};
            baseFlag = bitset(baseFlag,1);
        otherwise
            warning(['Unrecognized unit name: ' varargin{n} '. The corresponding value will go unused.'])
    end

end

% check if there is a power base input
if bitget(baseFlag, 7)
%   check if there is no frequency base  
    if ~bitget(baseFlag, 3)
        
%       check if there is voltage base input
        if bitget(baseFlag, 6) 
            i_base = s_base/v_base;
%           check if there is a current base input
            if bitget(baseFlag, 5)
                warning('Recaluclating Current base from Power and Voltage bases.')
            end
        
            z_base = v_base^2/s_base;
            if bitget(baseFlag, 4)
                warning('Recaluclating Impedance base from Power and Voltage bases.')
            end
%       if no voltage base then check if there is a current base
        elseif bitget(baseFlag, 5)
            v_base = s_base/i_base;
        
            z_base = s_base/i_base^2;
            if bitget(baseFlag, 4)
                warning('Recaluclating Impedance base from Power and Current bases.')
            end
%       if no voltage base and if no current base then check if there is an impedance base
        elseif bitget(baseFlag,4)
            v_base = sqrt(s_base*z_base);
            i_base = sqrt(s_base/z_base);
%   
        end
%   if there is a frequency base
    else
%       check if there is voltage base input
        if bitget(baseFlag, 6) 
            i_base = s_base/v_base;
%           check if there is a current base input
            if bitget(baseFlag, 5)
                warning('Recaluclating Current base from Power and Voltage bases.')
            end
        
            z_base = v_base^2/s_base;
            if bitget(baseFlag, 4)
                warning('Recaluclating Impedance base from Power and Voltage bases.')
            end
            
            c_base = 1/(z_base*2*pi*f_base);
            if bitget(baseFlag, 2)
                warning('Recaluclating Capacitance base.')
            end
            
            l_base = z_base/(2*pi*f_base);
            if bitget(baseFlag, 1)
                warning('Recaluclating Inductance base.')
            end
%       if no voltage base then check if there is a current base
        elseif bitget(baseFlag, 5)
            v_base = s_base/i_base;
        
            z_base = s_base/i_base^2;
            if bitget(baseFlag, 4)
                warning('Recaluclating Impedance base from Power and Current bases.')
            end
            
            c_base = 1/(z_base*2*pi*f_base);
            if bitget(baseFlag, 2)
                warning('Recaluclating Capacitance base.')
            end
            
            l_base = z_base/(2*pi*f_base);
            if bitget(baseFlag, 1)
                warning('Recaluclating Inductance base.')
            end
%       if no voltage base and if no current base then check if there is an impedance base
        elseif bitget(baseFlag,4)
            v_base = sqrt(s_base*z_base);
            i_base = sqrt(s_base/z_base);
            
%           c_base = 1/(z_base*2*pi*f_base);
            if bitget(baseFlag, 2)
                warning('Recaluclating Capacitance base.')
            end
            
            l_base = z_base/(2*pi*f_base);
            if bitget(baseFlag, 1)
                warning('Recaluclating Inductance base.')
            end
%       if no voltage base and if no current base and if no impedance base then check if there is a capacitance base
        elseif bitget(baseFlag,2)
            z_base = 1/(c_base*2*pi*f_base);
            v_base = sqrt(s_base*z_base);
            i_base = sqrt(s_base/z_base);
            
            l_base = z_base/(2*pi*f_base);
            if bitget(baseFlag, 1)
                warning('Recaluclating Inductance base.')
            end
%       if no voltage base and if no current base and if no impedance base and if no capacitance base then check if there is a iductance base
        elseif bitget(baseFlag,1)
            z_base = l_base*2*pi*f_base;
            v_base = sqrt(s_base*z_base);
            i_base = sqrt(s_base/z_base);
            c_base = 1/(z_base*2*pi*f_base);
        end

    end
    
% % no power base
% else        
end

% set the output based on the name of the input value and the bases
try
    switch lower(desired_unit)
        case powerNames
            unit_val = pu_val*s_base;
        case voltageNames
            unit_val = pu_val*v_base;
        case currentNames
            unit_val = pu_val*i_base;
        case impedanceNames
            unit_val = pu_val*z_base;
        case frequencyNames
            unit_val = pu_val*f_base;
        case capacitanceNames
            unit_val = pu_val*c_base;
        case inductanceNames
            unit_val = pu_val*l_base;
        otherwise
            warning(['Unrecognized unit name: ' desired_unit '. Could not calculate.'])
    end
catch
%     if the base was not calulated then there will be an error
   error('Could not calculate desired base and value from the inputs. Try adding additional base information.')
end
end


