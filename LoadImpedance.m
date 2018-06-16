function [ Z ] = LoadImpedance( P, Q, V )
%LOADIMPEDANCE 
%   This function takes the active (P) and reactive (Q) powers dissipated 
%   by a three phase load for a given line to line voltage (V) and 
%   calculates the equivalent per-phase impedance (Z).
% 
% The apparent power of a three phase system can be given as:
% S = 3*V_phase*conj(I) 
%   = 3*V_phase*conj(V_phase/Z)
%   = 3*V_phase*conj(V_phase)/conj(Z)
%   = 3*abs(V_phase)^2/conj(Z)
%   = abs(sqrt(3)*V_phase)^2/conj(Z)
%   = abs(V_line)^2/conj(Z)
% Solving for Z yeilds:
% Z = conj(abs(V_line)^2/S)

% Inputs:
%   P   Active 3 phase power in per-unit
%   Q   Reactive 3 phase power in per-unit
%   V   Line-to-line voltage magnitude in per-unit
% 
% Outputs;
%   Z   Equivalent phase impedance in per-unit

% Total apparent power 
S = P + 1i*Q;

% Caluclate Z using equation above
Z = conj(abs(V).^2./S);

end