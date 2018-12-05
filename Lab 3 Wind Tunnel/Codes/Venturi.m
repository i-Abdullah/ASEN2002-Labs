function [ Veloc Error ] = Venturi (P_atm, T_atm, AirP_Diff, sigma_P_atm, sigma_T_atm, sigma_Air_P_Diff,Rfluid,AreaRatio)
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
% and Venturi tube and ideal gas, the function will output the velocity and
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
%           AreaRatio : A2/A1 (Area of test / Area of the inital location)
%
%
% ----------------------- OUTPUT --------------------------------
%
%           Veloc : Velocity as matrix, in m/s or equivelent english units
%                   (if the inputs are in English units) 
%
%           Error : Uncertainties as matrix.
%
% ----------------------- IMPORTANT --------------------------------
% The only things treated as exact in error analysis are gas constant and
% area ratios.

%
% IMPORTANT: The order of the inputs is not arbitartry, it's because when
% MATLAB takes partial derivatives it re-arranges the inputs alphabetically
% and the easier sloution to have the uncertinity calculations to flip the
% order here rather than try to re-arrange the order of the derivatives.


%check if the user inputted un-equal arrays, because this code is meant to
%go and measure one data point or all a set of datas.

if length(P_atm) ~= length(T_atm)
    
    error('It seems that You atmospheric pressure and temp arrays are not equal, they must be equal, if you meant to have one of them constant repeat the value')
    
end

if length(P_atm) ~= length(AirP_Diff)
    
    error('It seems that You atmospheric pressure and pressure diffrence arrays are not equal, they must be equal, if you meant to have one of them constant repeat the value')
    
end

if length(T_atm) ~= length(AirP_Diff)
    
    error('It seems that You Temp and pressure diffrence arrays are not equal, they must be equal, if you meant to have one of them constant repeat the value')
    
end


%create the symbolic function


syms Pdiff P T
Velocity(Pdiff, P, T) = sqrt ( (2 * abs(Pdiff) * Rfluid * T) / ( (P) * (1-(AreaRatio)^2) ) ) ; 

%pre-define this matrix in case there's more than one Pressure or Temp
%measurments. 
Veloc = zeros(1,length(P_atm));
Error = zeros(1,length(P_atm));

%loop to find Velocity

for i = 1:length(P_atm)
Veloc(i) = double(Velocity(AirP_Diff(i),P_atm(i),T_atm(i)));
end


%% Error analysis: using partial derivatives method.

% Define Partial Derivatives.

    PartialDiffP = diff(Velocity,Pdiff);
    PartialPatm = diff(Velocity,P);
    PartialT= diff(Velocity,T);

    % loop to find error.
    
for i = 1:length(P_atm)

    Error(i) = (double(PartialDiffP(AirP_Diff(i),P_atm(i),T_atm(i)))*sigma_Air_P_Diff(i))^2 + ((double(PartialPatm(AirP_Diff(i),P_atm(i),T_atm(i))))*sigma_P_atm(i))^2 + ((double(PartialT(AirP_Diff(i),P_atm(i),T_atm(i)))*sigma_T_atm(i)))^2 ;

end

Error = sqrt(Error);



end
