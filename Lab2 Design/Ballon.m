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

AbsroptivitySun = 0.6;
qSun = 1370; % W / m^2 

% Albedo

qEarth = 237; % W / m^2
AbsroptivityEarth = 0.8;

% Emissivity of Ballon

EmissivityMaterial = 0.8;


%% Setup Equlibirum eqaution at the day and night

%-=-=-=-=-=-=-=-=-=-=-=( Find New T )=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

T_New_Day = ((( EmissivityMaterial * qEarth )/2 + (AbsroptivitySun*qSun)/2) / (4*EmissivityMaterial*StevBoltzConst))^(1/4);
T_New_Night = ((qEarth*AbsroptivityEarth*1/2)/(EmissivityMaterial*StevBoltzConst))^(1/4)


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

%-=-=-=-=-=-=-=-=-=-=-=( Find Density of Ballon )=-=-=-=-=-=-=-=-=-=-=-=

%New Density via Ideal Gas law ( 1 / Specific Volume).
NewDensityGas = ( ((PLoop+10)/1000) / (RGas*T_New_Night) );

% Solve the force Balance equation for the raduis of ballon
RaduisCuibed = MassLoad / ( (4*pi/3)  * ( rhoLoop - NewDensityGas - ( 3 * DensityMylar * ( (GagePressure * SafetyFactor) / (2*MueU) ) ) ) );
Raduis = RaduisCuibed^(1/3);
Thickness = ( (GagePressure*Raduis*SafetyFactor) / (2*MueU) );

VolumeBal = (4/3)*pi*RaduisCuibed;

DensityBallon = ( ( NewDensityGas*(4/3)*pi*(RaduisCuibed) ) + ( DensityMylar*4*pi*(Raduis)^2*Thickness) + (500) ) / ( ((4/3) * pi * RaduisCuibed) + (4*pi*(Raduis^2)*Thickness))


h = HuntHight(DensityBallon,0,80000);

while abs(rhoLoop-DensityBallon) > 1e-10
    
[ TLoop aLoop PLoop rhoLoop ] = atmoscoesa(h);

NewDensityGas = ( ((PLoop+10)/1000) / (RGas*T_New_Night) );
RaduisCuibed = MassLoad / ( (4*pi/3)  * ( rhoLoop - NewDensityGas - ( 3 * DensityMylar * ( (GagePressure * SafetyFactor) / (2*MueU) ) ) ) );
Raduis = RaduisCuibed^(1/3);
Thickness = ( (GagePressure*Raduis*SafetyFactor) / (2*MueU) );

VolumeBal = (4/3)*pi*RaduisCuibed;

DensityBallon = ( ( NewDensityGas*(4/3)*pi*(RaduisCuibed) ) + ( DensityMylar*4*pi*(Raduis)^2*Thickness) + (500) ) / ( ((4/3) * pi * RaduisCuibed) + (4*pi*(Raduis^2)*Thickness))

end

ActualHight = h;

%use ideal gas!

%% Overall Density: Night

[ TLoop aLoop PLoop rhoLoop ] = atmoscoesa(height);

%New Density via Ideal Gas law ( 1 / Specific Volume).
NewDensityGas_Night = ( ((PLoop+10)/1000) / (RGas*T_New_Night) );

% Solve the force Balance equation for the raduis of ballon
RaduisCuibed = MassLoad / ( (4*pi/3)  * ( rhoLoop - NewDensityGas_Night - ( 3 * DensityMylar * ( (GagePressure * SafetyFactor) / (2*MueU) ) ) ) );
Raduis = RaduisCuibed^(1/3);
Thickness = ( (GagePressure*Raduis*SafetyFactor) / (2*MueU) );

VolumeBal = (4/3)*pi*RaduisCuibed;

DensityBallon = ( ( NewDensityGas_Night*(4/3)*pi*(RaduisCuibed) ) + ( DensityMylar*4*pi*(Raduis)^2*Thickness) + (500) ) / ( ((4/3) * pi * RaduisCuibed) + (4*pi*(Raduis^2)*Thickness))

%%





%% Safet factor: Check if the stress will get up to the ultimat strength

k = 1; %safety factor
Pg = 10; %Gage Pressue in kPA

