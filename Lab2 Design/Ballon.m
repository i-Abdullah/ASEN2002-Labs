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

RaduisBallon = ;
SurfAreaBallon = 4 * pi * (r)^2 ;
VolumeBallon  = (4/3) * pi * (r)^3;


%% ?u = Qnet

%----------------------- ( Sources of Heat transfer ) ---------------------

% TEH SUNE:

Q_SUN_DAY = @(time) 1370 * time * (SurfAreaBallon/2);

% Albedo

Q_Albedo = @(time) 237 * time * (SurfAreaBallon/2) ;

% Emissivity of Ballon

EmissivityMaterial = 0.8;
StevBoltzConst = 5.670e-8;
Q_Emissivity = @(time,Temp) EmissivityMaterial * StevBoltzConst * SurfAreaBallon * Temp^4

% Constant Diffusuion?


%% Q = ?u ????? Q = Cv(T2 -T1) ????? Q/Cv + T1 = T2 ????? v = mRT2/P

%------------------ ( Find Total Q ) --------------------------------------
Time = ;
QTotal_Day = feval(Q_SUN_DAY,Time) - feval(Q_Emissivity,Time) + Q_Albedo*feval(Q_SUN_DAY,Time) ;
QTotal_Night = - feval(Q_SUN_Night,Time) - feval(Q_Emissivity,Time) - Q_Albedo*feval(Q_SUN_Night,Time); %Not sure about albedo at night

QTotat = [ QTotal_Day QTotal_Night ];

%------------------ ( Find New Volume ) --------------------------------------

% Find your initial info

height = 3500; %the initial hight

[ T a P rho ] = atmosisa(height);

RGas = 1; %Gas constant %CHEKKKKK!! ARE WE GOnna use funcitons at different temps?

T2 = T; %Temp at the original hight we start at 35 km.

P = 10 + P %Pressure inside our ballon taken to be 10 as gage, thus it's + Patm
CurrentVolume = massTotal*T2*R / P ;

%% Find the Density based on Volume

%% Find the new hight


%% Safet factor: Check if the stress will get up to the ultimat strength

k = 1; %safety factor
Pg = 10; %Gage Pressue in kPA

