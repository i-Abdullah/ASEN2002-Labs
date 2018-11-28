%% housekeeping

clear;
clc;
close all;

%% Read Data:

filename = 'sample_set_2018.csv';
Data = csvread(['Data/' filename ],1,0);

% naming

Patm = Data(:,1);
Tatm = Data(:,2);
Rohatm = Data(:,3);
airSpeed = Data(:,4);
pitotDynamicP = Data(:,5);
auxDynamicP = Data(:,6);

% Pressure ports: Scanivalve Pressure
SP1 = Data(:,7);
SP2 = Data(:,8);
SP3 = Data(:,9);
SP4 = Data(:,10);
SP5 = Data(:,11);
SP6 = Data(:,12);
SP7 = Data(:,13);
SP8 = Data(:,14);
SP9 = Data(:,15); % not connected
SP10 = Data(:,16);
SP11 = Data(:,17); % not connected
SP12 = Data(:,18);
SP13 = Data(:,19); % not connected
SP14 = Data(:,20);
SP15 = Data(:,21); % not connected
SP16 = Data(:,22);

angleOfAttack = Data(:,23);
stingNormForce = Data(:,24);
stingAxialForce = Data(:,25);
stingPitchingMoment = Data(:,26);

%% Extrapolate pressure:

%Explain:

%{

To do this we can extrapolate the pressure at the trailing edge (x = 3.5in)
based upon  pressure measured at the last two ports on both the upper and
lower wing surfaces (upper: p10, p8 and lower: p12, p14. These linear
extrapolations will give us two  different  estimates for the pressure
at the trailing edge, but we know that the pressures from the upper and
lower surfaces must converge at the trailing edge. Thus we should average
these  two pressure estimates to provide one estimate at the trailing  edge. 

%}

%Take Mean Values and linearly extrapolate.
%each angel of attack has 3 speeds

SP8_location = 2.1; %in m
SP10_location = 2.8; %in m
SP12_location = 2.8; %in m
SP14_location = 2.1; %in m

Ptrail = 0;

LocTrail = 3.5; %where the trail is in meters.

%each two conscutive ports on the same height will be in the same loop.
%so 10 and 8 together, and 12 and 14 together.

% SP 8 AND SP 10
for i=1:1:9
y = [ mean(SP8((i*20-19):i*20)); mean(SP10((i*20-19):i*20)) ];
t = [ SP8_location; SP10_location ];

%calculate the new pressure at the trailing edge.

Ptrail_UpperPorts(i) = Slope*(LocTrail) + Intercept;
Uncertinity_Ptrail_Upper(i) = sqrt( [ LocTrail 1] * Q * [ LocTrail ; 1 ] );

end

for i=1:1:9
y = [ mean(SP12((i*20-19):i*20)); mean(SP14((i*20-19):i*20)) ];
t = [ SP12_location; SP14_location ];
[ Slope Intercept SigY SigB SigM Q ] = LSM(t,y)

%calculate the new pressure at the trailing edge.

Ptrail_LowerPorts(i) = Slope*(LocTrail) + Intercept;
Uncertinity_Ptrail_Lower(i) = sqrt( [ LocTrail 1] * Q * [ LocTrail ; 1 ] );

end