%% HouseKeeping:
clear;
clc;
close all;

%% Define initial Values

MassLoad = 500; %kg
SafetyFactor = 2; %Safety Values
GagePressure = 10; % 10 pascals
MueU = 21.3e6 ; %Ultimate Tensile Strength Pa.
RGas = 2.0769; %Gas constant
StevBoltzConst = 5.670e-8; % Stevents Boltzman Constant
DensityMylar = 962; % kg/m^3 
HeliumSeaLevel = 0.164;

%-=-=-=-=-=-=-=-=-=-=-=( Heat transfer )=-=-=-=-=-=-=-=-=-=-=-=

%Small q is simply ?????

% TEH SUN:

AbsroptivitySun = 0.92;
qSun = 1370; % W / m^2 

% Albedo

qEarth = 237; % W / m^2
AbsroptivityEarth = 0.94;

% Emissivity of Ballon

EmissivityMaterial = 0.94;




%% 
%-=-=-=-=-=-=-=-=-=-=-=( Find raduis, and Gas Density )-=-=-=-=-=-=-=-=-=-=

% We will use a complex loop to first find the total new density, we will
% find the temp inside the ballon utilizing the heat equiliburim, use the
% new density and utilize a function we made to find the new hight, however
% the new density will be a summation of gas density and everything else
% density, so we will use the density of gas first to find the raduis,
% then use it to find the volume of the material, and get its density,
% then add them up.

%Get Initial hight, estimate Pressure there, then use it to find density

height = 35000;
[ TLoop aLoop PLoop rhoLoop ] = atmoscoesa(height);
NuetDensityGas = ( ((PLoop+10)/1000) / (RGas*TLoop) );



%-=-=-=-=-=-=-=-=-=-=-=( Find Density of Ballon )=-=-=-=-=-=-=-=-=-=-=-=

%New Density via Ideal Gas law ( 1 / Specific Volume).

% Solve the force Balance equation for the raduis of ballon

%@ 35 KM

RaduisCuibed = MassLoad / ( (4*pi/3)  * ( rhoLoop - NuetDensityGas - ( 3 * DensityMylar * ( (GagePressure * SafetyFactor) / (2*MueU) ) ) ) );
Raduis = RaduisCuibed^(1/3);
Thickness = ( (GagePressure*Raduis*SafetyFactor) / (2*MueU) );
VolumeShell = 4*pi*Thickness*(Raduis^2);
MassMaterial = VolumeShell * DensityMylar;
MassHeluim = NuetDensityGas * ((4/3) * pi * (RaduisCuibed));
TotalMass = MassHeluim+MassLoad+MassMaterial;

Volume35 = MassHeluim/NuetDensityGas;
VolumeGround = MassHeluim/HeliumSeaLevel;


%% Setup Equlibirum eqaution at the day and night

%-=-=-=-=-=-=-=-=-=-=-=( Find New T )=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

T_New_Day = ((( EmissivityMaterial * qEarth ) + (AbsroptivitySun*qSun)) / (4*EmissivityMaterial*StevBoltzConst))^(1/4);
T_New_Night = ((qEarth*AbsroptivityEarth)/(4*EmissivityMaterial*StevBoltzConst))^(1/4);


%% Day Calculations

% Find the new overall density, by first finding the new overall volume using
% the ideal gas law to find the volume and the pressure @ 35 km since it is
% not going to change much.

% get data @ 35 km.
height = 35000;
[ TLoop aLoop PLoop rhoLoop ] = atmoscoesa(height);

% Volume for ideal Gas
VolDay = (MassHeluim*RGas*T_New_Day)/((PLoop+10)/1000);

% Density from the constant toal mass
DensDay = (TotalMass)/VolDay;

%find the hight using a function we made.
h = HuntHight(DensDay,0,80000);

HightDay = h;


%% Night Calculations: starting from the hight we @ during day

height = HightDay;
[ TLoop aLoop PLoop rhoLoop ] = atmoscoesa(height);

% New Volume from ideal Gas Law.
VolNight = ( ( MassHeluim*RGas*(T_New_Night) ) / ((10+PLoop )/1000 ));

DensNight = TotalMass / VolNight ;

h = HuntHight(DensNight,0,80000);

HightNight = h;

%% Graphing
