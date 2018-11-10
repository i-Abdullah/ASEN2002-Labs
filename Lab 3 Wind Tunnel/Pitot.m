function [ Veloc Error ] = Pitot (P_atm, T_atm, AirP_Diff, sigma_P_atm, sigma_T_atm, sigma_Air_P_Diff,Rfluid)
% ASEN 2002, Lab 3: Wind Tunnel
%
%----------------------------------------------------------------
%
% Done By:
% 1- Jack Soltys
% 2- Foster Greer
% 3- Abdulla Al Ameri
% 4- Caelan Maitland
%
%----------------------------------------------------------------
%
% This function will utilize the airspeed equation using Bernoulli's equations
% and pito-static tube and ideal gas, the function will output the velocity and
% uncertinity calculations. It is making the assumption that the universal gas constant
% that comes from the ideal gas law is treated to be exact.
%
% ----------------------- INPUTS --------------------------------
%
%           P_atm : Atmospheric pressure, can be any size as long as it is
%           linear matrix, and everything else must be the same size. (Pa)
%           or equivelent english unit.
%
%           T_atm : Atmospheric Temperature, can be any size as long as it is
%           linear matrix, and everything else must be the same size. (K) or equivelent english unit.
%
%           AirP_Diff : Atmospheric Air Difference Pressure, can be any size as long as it is
%           linear matrix, and everything else must be the same size. (P) or equivelent english unit.
%
%           sigma_P_atm : Uncertainty in Atmospheric pressure, can be any size as long as it is
%           linear matrix, and everything else must be the same size.
%
%           sigma_T_atm : Uncertainty in Atmospheric Temperature, can be any size as long as it is
%           linear matrix, and everything else must be the same size.
%
%           sigma_Air_P_Diff : Uncertainty in Atmospheric Air Difference Pressure, can be any size as long as it is
%           linear matrix, and everything else must be the same size.
%
%           Rfluid : Universal constant of the fluid as gas, it must be in
%           (Pascal)(Pa), or equivelent english unit.
%
% ----------------------- OUTPUT --------------------------------
%
%           Veloc : Velocity as matrix, in m/s or equivelent english units
%                   (if the inputs are in English units) 
%
%           Error : Uncertainties as matrix.
%
% ----------------------- IMPORTANT --------------------------------
% The only thing treated as exact in error analysis is gas constant. 

%
% IMPORTANT: The order of the inputs is not arbitartry, it's because when
% MATLAB takes partial derivatives it re-arranges the inputs alphabetically
% and the easier sloution to have the uncertinity calculations to flip the
% order here rather than try to re-arrange the order of the derivatives.

syms Pdiff P T
Velocity(Pdiff, P, T) = sqrt ( (2 * Pdiff * Rfluid * T) / (P) ) ; 

%pre-define this matrix in case there's more than one Pressure or Temp
%measurments. 
Veloc = zeros(1,length(P_atm));
Error = zeros(1,length(P_atm));

%loop to find Velocity
for i = 1:length(P_atm)
Veloc(i) = double(Velocity(AirP_Diff,P_atm,T_atm));
end


% Error analysis: using partial derivatives method.

%Define Partial Derivatives.

    PartialDiffP = diff(Velocity,Pdiff);
    PartialPatm = diff(Velocity,P);
    PartialT= diff(Velocity,T);

    % loop to find error.
    
for i = 1:length(P_atm)

    Error(i) = (double(PartialDiffP(AirP_Diff,P_atm,T_atm))*sigma_Air_P_Diff)^2 + ((double(PartialPatm(AirP_Diff,P_atm,T_atm)))*sigma_P_atm)^2 + ((double(PartialT(AirP_Diff,P_atm,T_atm))*sigma_T_atm))^2 ;

end

Error = sqrt(Error);

end