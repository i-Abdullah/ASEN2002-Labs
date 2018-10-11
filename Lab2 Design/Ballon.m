%% HouseKeeping:
clear;
clc;
close all;

%% Measurmnets: Each measurment would be stored in an array, the first indicie is the measurment, the second is the uncertinty.


%---------------------- ( THE MASS ) ------------------------------------
%masses:
MassBallon = [ 9.80 0.005 ]; % in g
MassGas = [ 1.9304622 0.004 ]; % in g
MassLoad = [ 0.525 0.006 ]; % in g

%total masses:
TotalMass = MassBallon(1) + MassGas(1) + MassLoad(1) ;
TotalMass_Uncertintiy = MassBallon(2) + MassGas(2) + MassLoad(2);

% Fractional Masses:
FractMassGas = MassGas(1) / TotalMass ;
FractMassGas_Uncertinity = sqrt((MassGas(2)/MassGas(1))^2+ (TotalMass_Uncertintiy/TotalMass)^2) * FractMassGas;

FractMassLoad = MassLoad(1)/TotalMass;
FractMassLoad_Uncertinty = sqrt((MassLoad(2)/MassLoad(1))^2+ (TotalMass_Uncertintiy/TotalMass)^2) * FractMassLoad;

FractMassBallon = MassBallon(1)/TotalMass;
FractMassBallon_Uncertinity  = sqrt((MassBallon(2)/MassBallon(1))^2+ (TotalMass_Uncertintiy/TotalMass)^2) * FractMassBallon;

%--------------------------------------------------------------------------

%----------------- ( Volumes/Areas ) --------------------------------------------

RaduisBallon = 0.346075/2 ;
SurfAreaBallon = 4 * pi * (RaduisBallon)^2 ;
VolumeBallon  = (4/3) * pi * (RaduisBallon)^3;


%% ?u = Qnet

%----------------------- ( Sources of Heat transfer ) ---------------------

% TEH SUNE:

AbsroptivitySun = 0.6;
qSun = 1370 % W / m^2 

%Q_SUN_DAY = @(time) 1370 * time * (SurfAreaBallon/2);

% Albedo

qEarth = 237 % W / m^2

% Emissivity of Ballon

EmissivityMaterial = 0.8;
StevBoltzConst = 5.670e-8;

% Add up everything 


%% Setup Equlibirum eqaution at the new hight

%------------------ ( Find New T ) --------------------------------------

T_New = ((( EmissivityMaterial * qEarth )/2 + (AbsroptivitySun*qSun)/2) / (4*EmissivityMaterial*StevBoltzConst))^(1/4)

%------------------ ( Itteration ) --------------------------------------
% Find your initial info

RGas = 2.0769 ; %Gas constant
height = 35000; %the initial hight
[ TNew aNew PNew rhoNew ] = atmoscoesa(height);
NewDensity = ( (PNew/1000) / (RGas*T_Ne) ) )
h = WhatHight_roh(NewDensity);

while abs(rhoNew-NewDensity) > 1e-6
    
[ TNew aNew PNew rhoNew ] = atmoscoesa(h);


NewDensity = ( (PNew/1000) / RGas*T_New )
h = WhatHight_roh(NewDensity);

end

ActualHight = h;

%------------------ ( Find New Volume ) --------------------------------------

%% Find the Density based on Volume

%% Find the new hight


%% Safet factor: Check if the stress will get up to the ultimat strength

k = 1; %safety factor
Pg = 10; %Gage Pressue in kPA

