%% HouseKeeping:
clear;
clc;
close all;

%% call the data for Maylar and Poly before coating
BallonMaylar
BallonPolyE_NoCoating

%% Define initial Values

MassLoad = 500+100; %kg
SafetyFactor = 1.5; %Decided Based on Research
GagePressure = 10; % 10 pascals
MueU = 21.3e6 ; %Ultimate Tensile Strength Pa.
RGas = 2.0769; %Gas constant
StevBoltzConst = 5.670e-8; % Stevents Boltzman Constant
DensityPoly = 962; % kg/m^3 
HeliumSeaLevel = 0.164; %density of Heluim @ sea level

%-=-=-=-=-=-=-=-=-=-=-=( Heat transfer )=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

%Small q is amount of energy being transferred  

% TEH SUN:
AbsroptivitySun = 0.44; %unitless
qSun = 1370; % W / m^2 

% Albedo
qEarth = 237; % W / m^2
AbsroptivityEarth = 0.88; %unitless

% Emissivity of Ballon
EmissivityMaterial = 0.88; %unitless

%-=-=-=-=-=-=-=-=-=-=-=( Find New T )=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% Setup Equlibirum eqaution at the day and night to find the temp

T_New_Day = ((( EmissivityMaterial * qEarth ) + (AbsroptivitySun*qSun)) / (4*EmissivityMaterial*StevBoltzConst))^(1/4); % in Kelvin

T_New_Night = ((qEarth*AbsroptivityEarth)/(4*EmissivityMaterial*StevBoltzConst))^(1/4); % in Kelvin


%%

%-=-=-=-=-=-=-=-=-=-=-=( Find raduis, and Gas Density )-=-=-=-=-=-=-=-=-=-=


%Get Initial height, estimate Pressure there, then use it to find density

height = 35000;
[ TLoop aLoop PLoop rhoLoop ] = atmoscoesa(height);
NuetDensityGas = ( ((PLoop+10)/1000) / (RGas*TLoop) );


%-=-=-=-=-=-=-=-=-=-=-=( Initial State info )=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

RaduisCuibed = MassLoad / ( (4*pi/3)  * ( rhoLoop - NuetDensityGas - ( 3 * DensityPoly * ( (GagePressure * SafetyFactor) / (2*MueU) ) ) ) );
Raduis = RaduisCuibed^(1/3);
Thickness = ( (GagePressure*Raduis*SafetyFactor) / (2*MueU) );
VolumeShell = 4*pi*Thickness*(Raduis^2);
MassMaterial = VolumeShell * DensityPoly;
MassHeluim = NuetDensityGas * ((4/3) * pi * (RaduisCuibed));
TotalMass = MassHeluim+MassLoad+MassMaterial;
Volume35 = MassHeluim/NuetDensityGas;
VolumeGround = MassHeluim/HeliumSeaLevel;


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

DensNight = (TotalMass-100) / VolNight ;

h = HuntHight(DensNight,0,80000);

HightNight = h;

%%


%% Graphing

% The Following section will be for graphing. How some stuff will change
% with change in safety factor @ 35 km
i = 0;

for SF = 1:0.2:6
    i = i+1;
R(i) = MassLoad / ( (4*pi/3)  * ( rhoLoop - NuetDensityGas - ( 3 * DensityPoly * ( (GagePressure * SF) / (2*MueU) ) ) ) ); %RaduisCuibed
T(i) = ( (GagePressure*R(i)*SF) / (2*MueU) ); %ThicknessShell
VS(i) = 4*pi*T(i); %Volume Shell
MM(i) = VS(i) * DensityPoly; %Mass Material
end


figure(1)

subplot(2,2,1)
plot(R.^(1/3),1:0.2:6)
title('Radius vs Factor of Safety')
ylabel('FOS')
xlabel('Raduis (m)')
grid minor


subplot(2,2,2)
plot(T,1:0.2:6)
title('Thickness vs Factor of Safety')
ylabel('FOS')
xlabel('Thickness (m)')
grid minor


subplot(2,2,3)
plot(VS,1:0.2:6)
title('Shell Volume vs Factor of Safety')
ylabel('FOS')
xlabel('Volume (m^3)')
grid minor


subplot(2,2,4)
plot(MM,1:0.2:6)
title('Mass of the Material vs Factor of Safety')
ylabel('FOS')
xlabel('Mass (kg)')
grid minor



%%  

%establish the day and night values in a mtrix so we can comet plot them

figure(2)

x = 1:0.5:24;
y_coat = [ ones(1,23)*HightDay ones(1,24)*HightNight ];
y_Nocoat = [ ones(1,23)*HightDay_NoCoat ones(1,24)*HightNight_NoCoat ];
y_May = [ ones(1,23)*HightDayMay ones(1,24)*HightNightMay ];

scatter(x,y_coat,'^')
hold on
scatter(x,y_Nocoat,'*')
hold on
scatter(x,y_May,'o')

grid minor
ylim([33900 37000])
xlim([0 25])
box_x_1=[0 12 12 0];
box_y_1=[33900 33900 37000 37000];

box_x_2=[12 25 25 12];
box_y_2=[33900 33900 37000 37000];

patch(box_x_1,box_y_1,'green','FaceAlpha',0.08)
hold on
patch(box_x_2,box_y_2,'red','FaceAlpha',0.08)
hold on
ref = refline(0,35000)
ref.Color = 'b';
hold on
ref = refline(0,35800)
ref.Color = 'r';
hold on
ref = refline(0,34200)
ref.Color = 'r';

title('Altitude with different balloons')
xlabel('Time (hrs)')
ylabel('Altitude (m)')

legend('Poly with coating','Poly without coating','Mylar','Day cycle','Night cycle','Target Height','Upper bound','Lower bound')
