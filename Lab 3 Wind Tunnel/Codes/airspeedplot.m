%% 

%{

This code will do the plots for the airspeeds and voltages with the error
bars for all the groups.



%}


%%
clc;clear;close all

Voltageinitialpitot = [];
Velocitypitot = [];
Errorpitot = [];

Voltageinitialventuri = [];
Velocityventuri = [];
Errorventuri = [];

for i = 1:2:15
    if ismember(i,[1,5,9,13])
        [a1,b1,c1] = extractdata(11,i,1);
        Voltageinitialpitot = [Voltageinitialpitot;a1];
        Velocitypitot = [Velocitypitot,b1];
        Errorpitot = [Errorpitot,c1];
    else
        [a,b,c] = extractdata(11,i,2);
        Voltageinitialventuri = [Voltageinitialventuri;a];
        Velocityventuri = [Velocityventuri,b];
        Errorventuri = [Errorventuri,c];
    end
end

Voltagepitot = zeros(size(Voltageinitialpitot));
for j = 1:size(Voltageinitialpitot)
Voltagepitot(j) = str2double(Voltageinitialpitot{j});
end

Voltageventuri = zeros(size(Voltageinitialventuri));
for j = 1:size(Voltageinitialventuri)
Voltageventuri(j) = str2double(Voltageinitialventuri{j});
end
    
Velocitypitot = Velocitypitot';
Velocityventuri = Velocityventuri';

% slope and y-intercept for linear model
bmatpitot = [ones(size(Voltagepitot)),Voltagepitot];
mpitot = bmatpitot\Velocitypitot; 

bmatventuri = [ones(size(Voltageventuri)),Voltageventuri];
mventuri = bmatventuri\Velocityventuri; 

% plot
figure(1)
hold on
plot(Voltagepitot,mpitot(2)*Voltagepitot+mpitot(1))
errorbar(Voltagepitot,Velocitypitot,Errorpitot,'ro')
grid minor
title('Voltage vs. Velocity Measured by Airspeed Pressure Transducer')
xlabel('Voltage Supplied to Wind Tunnel (V)')
ylabel('Calculated Airspeed Velocity (m/s)')
legend('Linear Airspeed Model','Pitot-Static Probe Measurements','location','southeast')
xlim([0 10.5])
hold off

figure(2)
hold on
plot(Voltageventuri,mventuri(2)*Voltageventuri+mventuri(1))
errorbar(Voltageventuri,Velocityventuri,Errorventuri,'bo')
grid minor
title('Voltage vs. Velocity Measured by Airspeed Pressure Transducer')
xlabel('Voltage Supplied to Wind Tunnel (V)')
ylabel('Calculated Airspeed Velocity (m/s)')
legend('Linear Airspeed Model','Venturi Tube Measurements','location','southeast')
xlim([0 10.5])
hold off

function [Voltage,Velocityreturn,Errorreturn] = extractdata(SectionNum,GroupNum,exp)
%% define constants/hard code

RAir = 287.0 ; % pa / m^3 k
SigmaTemp = 0.25 ; % in k
SigmaatmPressure = (250-20)*10^3*(1.5/100); %from lab document
SigmaDiffPressure = 6894.76 * (1/100); %from lab document
AreaRatio = 1/9.5 ;

% input files
if GroupNum <= 9
    inputfileVV = ['VelocityVoltage_S0' num2str(SectionNum) '_G0' num2str(GroupNum) '.csv' ];
else
    inputfileVV = ['VelocityVoltage_S0' num2str(SectionNum) '_G' num2str(GroupNum) '.csv' ];
end

%% read the files/ call the reading fucntion : VV
[VV_Files,~] = SubsonicTunnel2(inputfileVV,'BoundaryLayer_S011_G01.csv',RAir);

%% Seperate and extract values
Data_VV = VV_Files(:,2);
VV_Patm = Data_VV{1,1}{:,1} ;
VV_Air_P_diff = Data_VV{2,1}{:,1};
VV_Aux_P_diff = Data_VV{3,1}{:,1};
VV_Temp = Data_VV{4,1}{:,1};
Voltage = Data_VV{6,1}{:,1};

[rows,~] = size(VV_Patm);

Patm_MeanValues_VV = zeros(rows/2,1);
Patm_stdValues_VV = zeros(rows/2,1);

Air_P_diff_MeanValues_VV = zeros(rows/2,1);
Air_P_diff_stdValues_VV = zeros(rows/2,1);

atmTemp_MeanValues_VV = zeros(rows/2,1);
atmTemp_stdValues_VV = zeros(rows/2,1);

Aux_P_diff_MeanValues_VV = zeros(rows/2,1);
Aux_P_diff_stdValues_VV = zeros(rows/2,1);

for i=1:((rows)/2)
    
    Patm_MeanValues_VV(i) = cell2mat(VV_Patm((2*i)-1,2));
    Patm_stdValues_VV(i) = cell2mat(VV_Patm((2*i),2));

    Air_P_diff_MeanValues_VV(i) = cell2mat(VV_Air_P_diff((2*i)-1,2));
    Air_P_diff_stdValues_VV(i) = cell2mat(VV_Air_P_diff((2*i),2));

    atmTemp_MeanValues_VV(i) = cell2mat(VV_Temp((2*i)-1,2));
    atmTemp_stdValues_VV(i) = cell2mat(VV_Temp((2*i),2));

    Aux_P_diff_MeanValues_VV(i) = cell2mat(VV_Aux_P_diff((2*i)-1,2));
    Aux_P_diff_stdValues_VV(i) = cell2mat(VV_Aux_P_diff((2*i),2));

end


%% asume constant uncertinity:
sigma_T_atm = ones(1,length(atmTemp_MeanValues_VV)) * SigmaTemp ;
sigma_Air_P_Diff = ones(1,length(Air_P_diff_MeanValues_VV)) * SigmaDiffPressure ;
sigma_P_atm = ones(1,length(Patm_MeanValues_VV)) * SigmaatmPressure ;

%% calculate Velocity
if exp==1
    [Velocityreturn,Errorreturn] = Pitot (Patm_MeanValues_VV, atmTemp_MeanValues_VV, Air_P_diff_MeanValues_VV, sigma_P_atm, sigma_T_atm,sigma_Air_P_Diff,RAir);
else
    [Velocityreturn,Errorreturn] = Venturi (Patm_MeanValues_VV, atmTemp_MeanValues_VV, Air_P_diff_MeanValues_VV, sigma_P_atm, sigma_T_atm, sigma_Air_P_Diff,RAir,AreaRatio);
end

end

function [ VV_Files,BL_Files ] = SubsonicTunnel2(inputfileVV,inputfileBL,RAir)
% ASEN 2002, Lab 3: Wind Tunnel, Fall 18.
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
% This function will read csv input files from CU Boulder ITLL subsonic
% wind tunnel and get the mean values for each voltage that the wind tunnel
% runs on, as well as standard deviation, and store them in cell arrays.
%
% ----------------------- INPUTS --------------------------------
%
%           inputfileVV : input csv file for volatage velocity experiment.
%
%           inputfileBB : input csv file for Boundary Layer experiment.
%
% ----------------------- OUTPUT --------------------------------
%
%           VV_Files : Velocity Voltage expirement reuslts, it will contain
%           a cell structure that's 4x2, Each one of the first columns will
%           tell you what data is in the next column (Second column). In the
%           second column you will have one of the 5 things: Information about Aux pressure difference,
%           Air pressure difference, Temperature, Atmospheric pressure, and Density
%           and  and  and Each one of the SECOND column will be another [ n x 2 ] array,
%           where n = number of different voltages the wind tunnel opereated on * 2 , and inside will be mean 
%           and standard deviation underneath each other in order.
%
%           BL_Files : since the Boundary Layer expirement opreats on one
%           voltage, but at different y-probe location (and fixed x
%           loaction bsed on the port of connection). This will also be
%           cell structure, the first column will tell you what is inside
%           the second column of cell arrays, and you'll have the
%           following: y-location of the probe, uncertainty in that
%           location, and finally the velocity of that location. Note that
%           the last one is the freestream velocity.


%{
% those are constants that are pre-defined, use them t
RAir = 287.0 ; % pa / m^3 k
SigmaTemp = 0.25 ; % in k, uncertainty in temperature readings. 
SigmaatmPressure = (250-20)*10^3*(1.5/100); %from lab document, uncertainty in atm pressure readings 
SigmaDiffPressure = 6894.76 * (1/100); %from lab document, uncertainty in Air pressure readings.
SigmaManometer = 0.1 ; % in inch, uncertainty from the readings of the manometer that's used for venturi tube expeirment. 
%}

%% read the data

%write the file names.
filename_VV = inputfileVV; %the vleocity voltage file name
filename_BL = inputfileBL; %the Boudnary layer file name


%check if file isn't csv then display error massege:

[ path1 name1 exten1 ] = fileparts(filename_VV) ;
[ path2 name2 exten2 ] = fileparts(filename_BL) ;

if strcmp(exten1, '.csv') == 0 || strcmp(exten2, '.csv') == 0
    
    error('At least on of your files are not csv files, please check them!');
    
end

%read

data_VV = csvread(['Data/' filename_VV],1,0);
data_BL = csvread(['Data/' filename_BL],1,0);

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

atm_P_BL = data_BL(:,1); %Atmospheric pressure in Pa.
atm_Temp_BL = data_BL(:,2); %Atmospheric Temperature in K.
air_diff_P_BL = data_BL(:,3); %Aire pressure Differntial in Pa.
Aux_diff_P_BL = data_BL(:,4); %Auxillery pressure Differntial in Pa.
Eld_x_BL = data_BL(:,5); %ELD Probe x axis location in mm.
Eld_y_BL = data_BL(:,6); %ELD Probe y axis location in mm.
Voltage_BL = data_BL(:,7); % Voltage data were recorded at (in Volts).

%% Dynamic naming

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% get how many different voltages we have in each file, and where they at.
%place holders

[ numvolt_VV location_VV ] = unique(Voltage_VV);
[ numvolt_BL location_BL ] = unique(Voltage_BL);

numvolt_VV = num2str(numvolt_VV);
numvolt_BL = num2str(numvolt_BL);

% the loop will check if there's digits in the voltage in replace them with
% underscore for the dynamic naming to work with any number.

for i =1:length(numvolt_VV)
    if contains(numvolt_VV(i),'.') ==1
        numvolt_VV(i) = strrep(a,'.','_')
    end
end

% repeat for BL file.

for i =1:length(numvolt_BL)
    if contains(numvolt_VV(i),'.') ==1
        numvolt_BL(i) = strrep(a,'.','_')
    end
end


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


%Dynamic naming

%pre-define the cells that will have the values

names_VV = {};
names_BL = {};

VV_atmP = {};
BL_atmP = {};

VV_Air_P_diff = {};
BL_Air_P_diff = {};

VV_Aux_P_diff = {};
BL_Aux_P_diff = {};

VV_Temp = {};
BL_Temp = {};

%the following loop is for Voltage Velocity file, and it'll get the mean
%values and STD at each voltage and name it accordingly. 
for i = 1:length(numvolt_VV)
    
    
 names_VV{i} = [num2str(numvolt_VV(i)) '_VV' ];
 
 % here the first row = mean, 2nd = std, and each column is different
 % voltage
 
 VV_atmP{(2*i)-1,1} = [ names_VV{i} '_mean_Patm'];
 VV_atmP{(2*i),1} = [ names_VV{i} '_std_Patm'];
 
 VV_Air_P_diff{(2*i)-1,1} = [ names_VV{i} '_mean_Air_P_diff'];
 VV_Air_P_diff{(2*i),1} = [ names_VV{i} '_std_Air_P_diff'];

 
 VV_Aux_P_diff{(2*i)-1,1} = [ names_VV{i} '_mean_Aux_P_diff'];
 VV_Aux_P_diff{(2*i),1} = [ names_VV{i} '_std_Aux_P_diff'];

 VV_Temp{(2*i)-1,1} = [ names_VV{i} '_mean_Temp'];
 VV_Temp{(2*i),1} = [ names_VV{i} '_std_Temp'];

 
 if i < length(numvolt_VV)
     %get mean values
     if location_VV(i) < location_VV(i+1)
 VV_Air_P_diff{(2*i)-1,2} = mean(air_diff_P_VV(location_VV(i):location_VV(i+1)-1)) ;
 VV_Temp{(2*i)-1,2} = mean(atm_Temp_VV( location_VV(i):location_VV(i+1)-1)) ;
 VV_Aux_P_diff{(2*i)-1,2} = mean(Aux_diff_P_VV( location_VV(i):location_VV(i+1)-1)) ;
 VV_atmP{(2*i)-1,2} = mean(atm_P_VV( location_VV(i):location_VV(i+1)-1)) ;
     else
 VV_Air_P_diff{(2*i)-1,2} = mean(air_diff_P_VV(location_VV(i+1)-1:location_VV(i))) ;
 VV_Temp{(2*i)-1,2} = mean(atm_Temp_VV( location_VV(i+1)-1:location_VV(i))) ;
 VV_Aux_P_diff{(2*i)-1,2} = mean(Aux_diff_P_VV( location_VV(i+1)-1:location_VV(i))) ;
 VV_atmP{(2*i)-1,2} = mean(atm_P_VV( location_VV(i+1)-1:location_VV(i))) ;
     end
    
    % get std values
    
 VV_Air_P_diff{(2*i),2} = std(air_diff_P_VV(location_VV(i):location_VV(i+1)-1)) ;
 VV_Temp{(2*i),2} = std(atm_Temp_VV( location_VV(i):location_VV(i+1)-1)) ;
 VV_Aux_P_diff{(2*i),2} = std(Aux_diff_P_VV( location_VV(i):location_VV(i+1)-1)) ;
 VV_atmP{(2*i),2} = std(atm_P_VV( location_VV(i):location_VV(i+1)-1)) ;


 else %condition for the last itteration
     
          %get mean values
 VV_Air_P_diff{(2*i)-1,2} = mean(air_diff_P_VV(location_VV(i):end)) ;
 VV_Temp{(2*i)-1,2} = mean(atm_Temp_VV(location_VV(i):end)); 
 VV_Aux_P_diff{(2*i)-1,2} = mean(Aux_diff_P_VV(location_VV(i):end)); 
 VV_atmP{(2*i)-1,2} = mean(atm_P_VV(location_VV(i):end));
    
    % get std values
    
 VV_Air_P_diff{(2*i),2} = std(air_diff_P_VV((location_VV(i):end))); 
 VV_Temp{(2*i),2} = std(atm_Temp_VV((location_VV(i):end)));
 VV_Aux_P_diff{(2*i),2} = std(Aux_diff_P_VV((location_VV(i):end)));
 VV_atmP{(2*i),2} = std(atm_P_VV((location_VV(i):end)));
 

 end
 

end
 


%% density values

%density for VV and BL (Velocity Voltage and Boundary layer, each will be done in
%a seperate loop.


% we will store all the matrices in a giant matrix, each row represents the density values at that
% voltage.

VV_Density = zeros(length(Voltage_VV),1);
BL_Density = zeros(length(Voltage_BL),1);


% we will run a loop to compute the densities and place them

% each input file will have its own loop


% the vleocity voltage loop
for i=1:length(numvolt_VV)
    
    %condition for the last itteration
    if i==length(numvolt_VV)
        
        for j = location_VV(i):length(Voltage_VV)
        VV_Density(j) = atm_P_VV(j)/( atm_Temp_VV(j) * RAir) ;
        end
        
    else
        
    for j = location_VV(i):location_VV(i+1)
        
        VV_Density(j) = atm_P_VV(j)/( atm_Temp_VV(j) * RAir) ;
    end
    
    end
    
end


% the boundary layer loop

for i=1:length(numvolt_BL)
    
    %condition for the last itteration
    if i==length(numvolt_BL)
        
        for j = location_BL(i):length(Voltage_BL)
        BL_Density(j) = atm_P_BL(j)/( atm_Temp_BL(j) * RAir) ;
        end
        
    else
        
    for j = location_BL(i):location_BL(i+1)
        
        BL_Density(j) = atm_P_BL(j)/( atm_Temp_BL(j) * RAir) ;
    end
    
    end
    
end


%% get the BL information

%we will need to know where (what row of the file) the data switches
%location, and in total there's 12 locations where the data will switch,
%where the y-probe will move, so we will look @ column #6.

%{
ChangeIndi = {} ; %change index array
j = 1;
for i = 2:length(Eld_y_BL);
    
    if abs(abs(Eld_y_BL(i))-abs(Eld_y_BL(i-1))) > 0.01
    ChangeIndi{j} = i;
    j = j+1;
    
    end
end
ChangeIndi = [ 1 , cell2mat(ChangeIndi) ];
BL_YprobeLocation = {}; %location of the y-probe, note that the last on is free-stream condition
BL_SigmaYprobeLocation = {}; %uncertainty in that location via std note that the last on is free-stream condition
BL_Velocity = {} ; %Velocity at each y-probe location, and note that the last on is free-stream condition
for i = 1:12
    
    if i==12
        
BL_YprobeLocation{i,1} = mean(Eld_y_BL(ChangeIndi(i):end));
BL_SigmaYprobeLocation{i,1} = std(Eld_y_BL((ChangeIndi(i):end)));
BL_Velocity{i,1} = sqrt(( 2 * mean(Aux_diff_P_BL((ChangeIndi(i):end)) * RAir * mean(atm_Temp_BL( ( ChangeIndi(i):end) )) / (mean(atm_P_BL(ChangeIndi(i):end))))));
    else
BL_YprobeLocation{i,1} = mean(Eld_y_BL(ChangeIndi(i):ChangeIndi(i+1)));
BL_SigmaYprobeLocation{i,1} = std(Eld_y_BL((ChangeIndi(i):ChangeIndi(i+1))));
BL_Velocity{i,1} = sqrt(( 2 * mean(Aux_diff_P_BL((ChangeIndi(i):ChangeIndi(i+1))) * RAir * mean(atm_Temp_BL( ( ChangeIndi(i):ChangeIndi(i+1)) )) / ( mean(atm_P_BL(ChangeIndi(i):ChangeIndi(i+1)))))));
    end
end
%}

%another way that is not atonomus :

% the data are parsed at every 500, each one of these will be [ 3 x 12 ]

% 12 : 12 different voltages and the center is the 12th.

% 3 : the first row is the mean value and the second is std, third is
% velocity at that location.


BL_YprobeLocation = {}; %location of the y-probe, note that the last on is free-stream condition
BL_SigmaYprobeLocation = {}; %uncertainty in that location via std note that the last on is free-stream condition
BL_Velocity = {} ; %Velocity at each y-probe location, and note that the last on is free-stream condition

for i = 1:12
    
BL_YprobeLocation{i,1} = mean(Eld_y_BL(i*500-499:i*500));
BL_SigmaYprobeLocation{i,1} = std(Eld_y_BL(i*500-499:i*500));
BL_Velocity{i,1} = sqrt(( 2 * mean(Aux_diff_P_BL(i*500-499:i*500)) * RAir * mean(atm_Temp_BL(i*500-499:i*500)) ) / ( mean(atm_P_BL((i*500-499:i*500)) )));

end



%% store VV values

VV_Files = {};

VV_Files{1,1} = { 'VV_atmP' };
VV_Files{2,1} = { 'VV_Air_P_diff' };
VV_Files{3,1} = { 'VV_Aux_P_diff' };
VV_Files{4,1} = { 'VV_Temp' };
VV_Files{5,1} = { 'VV_Density' };
VV_Files{6,1} = { 'Voltages' };

VV_Files{1,2} = { VV_atmP };
VV_Files{2,2} = { VV_Air_P_diff };
VV_Files{3,2} = { VV_Aux_P_diff };
VV_Files{4,2} = { VV_Temp };
VV_Files{5,2} = { VV_Density };
VV_Files{6,2} = { cellstr(numvolt_VV) };

%% store BL values

BL_Files = {};

BL_Files{1,1} = { 'BL_YprobeLocation' };
BL_Files{2,1} = { 'BL_SigmaYprobeLocation' };
BL_Files{3,1} = { 'BL_Velocity' };

BL_Files{1,2} = { BL_YprobeLocation };
BL_Files{2,2} = { BL_SigmaYprobeLocation };
BL_Files{3,2} = { BL_Velocity };

end