%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this MATLAB code is to read data files collected in the
%%% wind tunnel. This MATLAB script is able to quantify as well as analyze
%%% things like the airspeed in a the test section using different methods,
%%% the errors that occur when collecting data, and the thickeness of the
%%% boundary layer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Date created: Dec/24/2018
%%% Last modefied: Nov/12/2018
%%% Authors: Stephen Chamot & Abdullah Almugairin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Assumption: Incompressible flow, Ideal Gas laws apply
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping

clear all;
close all;
clc;
%% Read data

% Reading velocity voltage data 
dataRawV = csvread('VelocityVoltage_S011_G13.csv',1,0);

% Create matrices that will be used later 
matrix1 = [];
matrix2 = [];
matrix3 = [];
matrix4 = [];
matrix5 = [];

% Reading manometer data 

ManoData = csvread('Manometer Readings - Sheet1.csv',1,0);
ManoRdngDiff = ManoData(1:end,4);
ManoDiffPressue = ManoRdngDiff.*248.84;

ManoDiffPressue1 = ones(500,1) * ManoDiffPressue(1);
ManoDiffPressue2 = ones(500,1) * ManoDiffPressue(2);
ManoDiffPressue3 = ones(500,1) * ManoDiffPressue(3);
ManoDiffPressue4 = ones(500,1) * ManoDiffPressue(4);
ManoDiffPressue5 = ones(500,1) * ManoDiffPressue(5);

ManoDiffPressueTotal = [ManoDiffPressue1;ManoDiffPressue2;ManoDiffPressue3;ManoDiffPressue4;ManoDiffPressue5];
%{
% Manometer readings. Multilply by 248.84 to convert from inH2o to Pa
ManoDiffPressue1 = ones(500,1) *  0.05.*248.84;
ManoDiffPressue2 = ones(500,1) * 0.43.*248.84;
ManoDiffPressue3 = ones(500,1) * 1.5.*248.84;
ManoDiffPressue4 = ones(500,1) * 2.9.*248.84;
ManoDiffPressue5 = ones(500,1) * 4.9.*248.84;
%}
%% Find voltages used from Velocity Voltage data file

% Get the rows and cpolumns
[rows, columns] = size(dataRawV);

% Assign constants
i = 1;
j = 1;

% Repeat the process 5 times for each of the 5 voltages
volt1 = dataRawV(i, 7); 

% While loop to find the first voltage
while dataRawV(i, 7) == volt1
     i = i + 1;
end
% Save all the data in a matrix when the first voltage was used
matrix1 = dataRawV(j:i-1,:);
j = i;

% While loop to find the second voltage
volt2 = dataRawV(i, 7); 
while dataRawV(i, 7) == volt2
     i = i + 1;
end
% Save all the data in a matrix when the second voltage was used
matrix2 = dataRawV(j:i-1,:);
j = i;

% While loop to find the third voltage
volt3 = dataRawV(i, 7); 
while dataRawV(i, 7) == volt3
     i = i + 1;
end
% Save all the data in a matrix when the third voltage was used
matrix3 = dataRawV(j:i-1,:);
j = i;

% While loop to find the fourth voltage
volt4 = dataRawV(i, 7); 
while dataRawV(i, 7) == volt4
     i = i + 1;
end
% Save all the data in a matrix when the fourth voltage was used
matrix4 = dataRawV(j:i-1,:);
j = i;

% While loop to find the fifth voltage
volt5 = dataRawV(i, 7); 
while dataRawV(i, 7) == volt5
     i = i + 1;
     if rows < i
         break;
     end
end
% Save all the data in a matrix when the fifth voltage was used
matrix5 = dataRawV(j:i-1,:);
j = i;

% Put voltages in one column 
Volts = [volt1;volt2;volt3;volt4;volt5];

%% Compute velocity from pitot probe and static ring

% Constants
R = 287; % gas constant of air in J kg^-1 K^-1
A1A2ratio = 1/9.5; % Area ratio

% Airspeed differential pressure, temperature and pressure
%for the first voltage
diffPressure1 = matrix1(:,3);
temp1 = matrix1(:,2);
patm1  = matrix1(:,1);

% Calculate the velocity for the first voltage for pitot probe and static ring
V1 = sqrt((2 .* diffPressure1 .* R .* temp1 ) ./ (patm1 .* (1 - (A1A2ratio)^2)));
VentV1 = sqrt((2 .* ManoDiffPressue1 .* R .* temp1 ) ./ (patm1 .* (1 - (A1A2ratio)^2)));


% Airspeed differential pressure, temperature and pressure
%for the second voltage
diffPressure2 = matrix2(:,3); 
temp2 = matrix2(:,2);
patm2  = matrix2(:,1);

% Calculate the velocity for the second voltage for pitot probe and static ring
V2 = sqrt((2 .* diffPressure2 .* R .* temp2 ) ./ (patm2 .* (1 - (A1A2ratio)^2)));
VentV2 = sqrt((2 .* ManoDiffPressue2 .* R .* temp2 ) ./ (patm2 .* (1 - (A1A2ratio)^2)));


% Airspeed differential pressure, temperature and pressure
%for the third voltage
diffPressure3 = matrix3(:,3); %%Airspeed diff pressure
temp3 = matrix3(:,2);
patm3  = matrix3(:,1);

% Calculate the velocity for the third voltage for pitot probe and static ring
V3 = sqrt((2 .* diffPressure3 .* R .* temp3 ) ./ (patm3 .* (1 - (A1A2ratio)^2)));
VentV3 = sqrt((2 .* ManoDiffPressue3 .* R .* temp3 ) ./ (patm3 .* (1 - (A1A2ratio)^2)));

% Airspeed differential pressure, temperature and pressure
%for the fourth voltage
diffPressure4 = matrix4(:,3); %%Airspeed diff pressure
temp4 = matrix4(:,2);
patm4  = matrix4(:,1);

% Calculate the velocity for the fourth voltage for pitot probe and static ring
V4 = sqrt((2 .* diffPressure4 .* R .* temp4 ) ./ (patm4 .* (1 - (A1A2ratio)^2)));
VentV4 = sqrt((2 .* ManoDiffPressue4 .* R .* temp4 ) ./ (patm4 .* (1 - (A1A2ratio)^2)));

% Airspeed differential pressure, temperature and pressure
%for the fifth voltage
diffPressure5 = matrix5(:,3); %%Airspeed diff pressure
temp5 = matrix5(:,2);
patm5  = matrix5(:,1);

% Calculate the velocity for the fifth voltage for pitot probe and static ring
V5 = sqrt((2 .* diffPressure5 .* R .* temp5 ) ./ (patm5 .* (1 - (A1A2ratio)^2)));
VentV5 = sqrt((2 .* ManoDiffPressue5 .* R .* temp5 ) ./ (patm5 .* (1 - (A1A2ratio)^2)));

% Concatenate velocities 
Vpitot = [V1;V2;V3;V4;V5];
VenturiVelocity = [VentV1;VentV2;VentV3;VentV4;VentV5];


%% Compute the uncertainty in the two velocities

% Call the function that calculates the uncertainties in the velocities
% Get the air differential, P and T atm
Airdiff = dataRawV(1:end,3);
Tatm = dataRawV(1:end,2);
Patm = dataRawV(1:end,1);
x = size(Patm);
% Uncrtainties in tools
sigmaDP = 6894.76*(1/100);
sigmaP = (250-20)*((1.5)/100);
sigmaT = .25;

sigmaD_P = ones(1,length(ManoDiffPressueTotal))*sigmaDP;
sigma_P = ones(1,length(ManoDiffPressueTotal))*sigmaP;
sigma_T = ones(1,length(ManoDiffPressueTotal))*sigmaT;



[sigmaPitot] = SigmaV(R,ManoDiffPressueTotal,Tatm,Patm,x,sigmaD_P,sigma_P,sigma_T);
[sigmaVent] = SigmaVelVent(R,A1A2ratio,Airdiff,Tatm,Patm,x,sigmaD_P,sigma_P,sigma_T);







%% Plots for the velocities, error bars and pressures

% Plot the velocities against each other and label
figure (1)
time = 1:1:2500;
plot(time,Vpitot,'r.',time,VenturiVelocity,'b.')
xlabel('Time (s)')
ylabel('Airspeed (m/s)')
legend('Water Manometer - Pitot-Static Probe','Pressure Transducer - Static Ring','location','northwest')
title('Group 13 Airspeed vs. Time Using Different Measurement Devices')

% Plot and label the figure for the uncertainties in the Pitot tube
% measurements
VpitotMean = [mean(V1);mean(V2);mean(V3);mean(V4);mean(V5)];

figure (2)
subplot(1,2,1)
plot(Volts,VpitotMean,'o')
hold on
title('Voltage vs Average Airspeed (Using Pitot Tube) With a Best Fit')
xlabel('Voltage, Volts')
ylabel('Average Airspeed (m/s)')
legend('Average Airspeed','location','northwest')
xlim([0 12])
hold on
% Best fit line
coeff = polyfit(Volts,VpitotMean,1); %coefficients 
bfLine = coeff(1).*Volts + coeff(2); %data points
plot(Volts,bfLine) %plot voltages vs best fit line
legend('Average Airspeed','Best-Fit Line','location','northwest')
 
% Plot and label the figure for the uncertainties in the Venturi
% measurements
VVentmean = [mean(VentV1);mean(VentV2);mean(VentV3);mean(VentV4);mean(VentV5)];
subplot(1,2,2)
plot(Volts,VVentmean,'o')
title('Voltage vs Average Airspeed (Using Static Ring) With a Best Fit')
xlabel('Voltage (V)')
ylabel('Average Airspeed (m/s)')
legend('Average Airspeed','location','northwest')
xlim([0 12])
hold on
% Best fit line
coeff = polyfit(Volts,VVentmean,1); %coefficients 
bfLine = coeff(1).*Volts + coeff(2); %data points
plot(Volts,bfLine) %plot voltages vs best fit line
legend('Average Airspeed','Best-Fit Line','location','northwest')

% Error in the pitot probe velocity measurments 
figure (3)
subplot(2,1,1)
err = mean(sigmaVPitot);
err = err';
errorbar(VVentmean,Volts,err,'.','horizontal')
hold on
title('Error bars in velocity measurments of Venturi Tube')
xlabel('Airspeed (m/s)')
ylabel('Voltage (V)')
legend('Error Bars','location','northwest')

% Error in the manometer velocity measurments 
subplot(2,1,2)
err2 = mean(sigmaVVent);
err2 = err2';
errorbar(VpitotMean,Volts,err2,'.','horizontal')
title('Error bars in velocity measurments of Pitot-Static Probe')
xlabel('Airspeed (m/s)')
ylabel('Voltage (V)')
legend('Error Bars','location','northwest')

% Plot the pressures of the manometer and the pitot static
figure(4)
ManoDiffPressueTotal = [ManoDiffPressue1;ManoDiffPressue2;ManoDiffPressue3;ManoDiffPressue4;ManoDiffPressue5];
diffPressureTotal = [diffPressure1;diffPressure2;diffPressure3;diffPressure4;diffPressure5];
time = 1:1:2500;
plot(time,ManoDiffPressueTotal,'r.',time,diffPressureTotal,'b.')
title('Differential pressure from U-tube manometer and pressure tranducer')
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('Pressure from U-tube Manometer','Pressure from pressure pranducer','location','northwest')


%% Boundry layer analysis 

% Get and saperate data
BLData = csvread('BoundaryLayer_S011_G13.csv',1,0);
pressureB = mean(BLData(:,1));
tempB = mean(BLData(:,2));
auxDiffPressure = BLData(:,4);
ELDy = BLData(:,6).*10^(-3); % Convert to meters
BLvolt = BLData(1:end,7);

% Calculate the centerline airspeed 
BLvel = sqrt(2.*auxDiffPressure.*R.*tempB./pressureB);

%centerline velocity average
CnterlineVel = mean(BLvel(5502:end));

% Compute the anitial area and convert it to meters squared
TestSecArea1 = 12*12;
TestSecArea1 = TestSecArea1.*0.00064516; 

% Get coefficients and find the area inside the test section
BLVel = coeff(1)*5 + coeff(2);
TestSecArea2 = TestSecArea1*BLVel/CnterlineVel;

% While loop to calculate the edge of the bounry layer where the
% free-stream velocity reaches 95%
for i=1:length(BLvel)
    if BLvel(i) >= BLvel*0.95
        break;
    end
end

% Thickness of the boundry layer in meters squared
thickness = ELDy(i);



%% Compute the theoretical the boundary layer for a turbulent and a laminar flow

%air density
pair = 1.225; %kg/m^3
%length of test section
length = .6096; %meters
%viscocity of air at standard atmospheric temperature
viscocity = 1.88*10^-5;
%kinematic viscocity
KinViscocity = viscocity/pair;

%estimate the boundary layer for turbulent and laminar flow
dturbulent = 0.37*length^.8*(KinViscocity/BLVel)^.2; 
dlaminar = 5.2*sqrt(KinViscocity*length/BLVel);

%expected centerline velocity change
expV2 = TestSecArea1*BLVel/(.3048-dlaminar)^2;

%A2 found using calculated thickness from experimental data
expA2 = (sqrt(TestSecArea1)-thickness)^2;

%% Output Values

sigma_Pitot = {sigmaPitot(1,1:500);sigmaPitot(1,501:1000);sigmaPitot(1,1001:1500)...
    ;sigmaPitot(1,1501:2000);sigmaPitot(1,2001:2500)};


sigmaVent = {sigmaVent(1,1:500);sigmaVent(1,501:1000);sigmaVent(1,1001:1500)...
    ;sigmaVent(1,1501:2000);sigmaVent(1,2001:2500)};

Outputs = table(VpitotMean,sigma_Pitot,VpitotMean,sigmaVent)

fprintf('The theoretical thickness for a laminar flow is %f meters',dturbulent)
fprintf('The theoretical thickness for a laminar flow is %f meters',dlaminar)
