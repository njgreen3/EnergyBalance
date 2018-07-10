% This script tests the VSI_semiconductor_loss.m function


Pload = 5000;    % W
Qload = 2000;    % VAR
Vac = 480;   % V
switchFreq = 50e3;   % Hz
Vthresh = 4; % V
Vblock = 5;   % V
Rcond = 1e-3; % ohm
turnOnTime = 3e-6;   % s
turnOffTime = 3e-6;  % s
duty = .5;   
modAmp = .8;
swTheta1 = 0; %?
swTheta2 = .1*pi; %?  
condTheta1 = (0:.01:1)*pi;   %?
% condTheta2 = 1*pi;   %?
condTheta2 = condTheta1 + .1*pi;

Pcond = zeros(size(condTheta2));
Psw = zeros(size(condTheta2));

for n=1:length(condTheta2)
[ Pcond(n), Psw(n) ] = VSI_semiconductor_loss(Pload, Qload, Vac, switchFreq, ...
                        Vthresh, Vblock, Rcond, turnOnTime, turnOffTime, ...
                        duty, modAmp, swTheta1, swTheta2, condTheta1(n), condTheta2(n) );
                 
end


plot(condTheta1,Pcond)
xlabel('Conduction angle start (rad)')
ylabel('Conduction loss (W)')
