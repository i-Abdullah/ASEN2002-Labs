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

%% read the data

%write the file names.

filename_VV = 'VelocityVoltage_S011_G01.csv'; %the vleocity voltage file name
filename_BL = 'BoundaryLayer_S011_G01.csv'; %the Boudnary layer file name


%read

data_VV = csvread(filename_VV,1,0);
data_BL = csvread(filename_BL,1,0);

% divide, and conquer.

%from Velocity Voltage data
atm_P_VV = data_VV(:,1); %Atmospheric pressure in Pa.
atm_Temp_VV = data_VV(:,2); %Atmospheric Temperature in K.
air_diff_P_VV = data_VV(:,3); %Aire pressure Differntial in Pa.
Aux_diff_P_VV = data_VV(:,4); %Auxillery pressure Differntial in Pa.
Eld_x_VV = data_VV(:,5); %ELD Probe x axis location in mm.
Eld_y_VV = data_VV(:,6); %ELD Probe y axis location in mm.
Voltage_VV = data_VV(:,7); % Voltage data were recorded at (in Volts).

%frm Boundary layer

atm_P_BB = data_BL(:,1); %Atmospheric pressure in Pa.
atm_Temp_BB = data_BL(:,2); %Atmospheric Temperature in K.
air_diff_P_BB = data_BL(:,3); %Aire pressure Differntial in Pa.
Aux_diff_P_BB = data_BL(:,4); %Auxillery pressure Differntial in Pa.
Eld_x_BB = data_BL(:,5); %ELD Probe x axis location in mm.
Eld_y_BB = data_BL(:,6); %ELD Probe y axis location in mm.
Voltage_BB = data_BL(:,7); % Voltage data were recorded at (in Volts).

%%


