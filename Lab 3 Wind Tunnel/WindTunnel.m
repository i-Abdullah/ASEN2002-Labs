%% info

%{

Thie script will extract and analyze data obtained from a wind tunnel lab,
part of ASEN 2002: Lab 3, CU Bouldr, Fall 18.


Done by

1- Jack Soltys
2- Abdulla Al Ameri
3- Greer Foster
4- Caelan Maitland


%}



%% housekeeping

clear;
clc;
close all;

%% define constants/ hard code

% those are constants, although 
RAir = 287.0 ; % pa / m^3 k
SigmaTemp = 0.25 ; % in k
SigmaatmPressure = (250-20)*10^3*(1.5/100); %from lab document
SigmaDiffPressure = 6894.76 * (1/100); %from lab document
SigmaManometer = 0.1 ; % in inch

inputfileVV = 'VelocityVoltage_S011_G01.csv';
inputfileBL = 'BoundaryLayer_S011_G01.csv';

%% read the files/ call the reading fucntion

[ VV_Files BL_Files ] = SubsonicTunnel(inputfileVV,inputfileBL,RAir);

%extract values:
Data_VV = VV_Files(:,2);;
VV_Patm = Data_VV{1,1}{:,1} ;
VV_Air_P_diff = Data_VV{2,1}{:,1};
VV_Aux_P_diff = Data_VV{3,1}{:,1};
VV_Temp = Data_VV{4,1}{:,1};
VV_Density = Data_VV{5,1}{:,1};
VV_Voltages = Data_VV{6,1}{:,1};


% put all the mean values together, all the std values together.

[ rows columns ] = size(VV_Patm);

Patm_MeanValues = zeros(rows,1);
Patm_stdValues = zeros(rows,1);

for i=1:length(rows)/2
    
    Patm_MeanValues(i) = double(VV_Patm(i,2));
    Patm_stdValues(i) = double(VV_Patm(i+1,2));

end

%repeat for each file

[ rows columns ] = size(VV_Air_P_diff);

Air_P_diff_MeanValues = zeros(rows,1);
Air_P_diff_stdValues = zeros(rows,1);

for i=1:length(rows)/2
    
    Air_P_diff_MeanValues(i) = double(VV_Patm(i,2));
    Air_P_diff_stdValues(i) = double(VV_Patm(i+1,2));

end


%}

%% printout the results:
%{
Voltage = VV_Files{6,2}{:,1};
Error_Vento = { error_Vento_1 ; error_Vento_3 ; error_Vento_5 ; error_Vento_7 ; error_Vento_9 };
Veloc_Pitot = { Velocity_VV_1_Pito ; Velocity_VV_3_Pito ; Velocity_VV_5_Pito ; Velocity_VV_7_Pito ; Velocity_VV_9_Pito};
Error_BL = {'N/A';'N/A';error_BL_5;'N/A';'N/A'};
Error_Pitot = { error_pitot_1 ; error_pitot_3 ; error_pitot_5 ; error_pitot_7 ; error_pitot_9 };
Veloc_Venturi = { Velocity_VV_1_Vento ; Velocity_VV_3_Vento ; Velocity_VV_5_Vento ; Velocity_VV_7_Vento ; Velocity_VV_9_Vento};
Veloc_BL = { 'N/A' ; 'N/A' ; Velocity_BL_5_Pito ; 'N/A' ; 'N/A' };

Results = table(Voltage,Veloc_Pitot,Error_Pitot,Veloc_Venturi,Error_Vento,Veloc_BL,Error_BL)

%% Boundary Layer:

%% Plots

figure(1);
plot(1,Velocity_VV_1_Pito,'*')
hold on
errorbar(1,Velocity_VV_1_Pito,error_pitot_1);
hold on
plot(1.05,Velocity_VV_1_Vento,'*')
hold on
errorbar(1.05,Velocity_VV_1_Vento,error_Vento_1);
xlim([0 2]);
legend('','Ptio-stat, V=1','','Venturi, V=1');
grid minor

%}