function [ Results ] = SubsonicTunnel(inputfileVV,inputfileBB,ManoReadings)
%
%
%
%
%
%
%
%
%
%
%
%
%
%% define constants/ hard code

% those are constants that are pre-defined, use them t
RAir = 287.0 ; % pa / m^3 k


SigmaTemp = 0.25 ; % in k, uncertainty in temperature readings. 
SigmaatmPressure = (250-20)*10^3*(1.5/100); %from lab document, uncertainty in atm pressure readings 
SigmaDiffPressure = 6894.76 * (1/100); %from lab document, uncertainty in Air pressure readings.
SigmaManometer = 0.1 ; % in inch, uncertainty from the readings of the manometer that's used for venturi tube expeirment. 
%% read the data

%write the file names.
filename_VV = inputfileVV; %the vleocity voltage file name
filename_BL = inputfileBB; %the Boudnary layer file name

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

atm_P_BL = data_BL(:,1); %Atmospheric pressure in Pa.
atm_Temp_BL = data_BL(:,2); %Atmospheric Temperature in K.
air_diff_P_BL = data_BL(:,3); %Aire pressure Differntial in Pa.
Aux_diff_P_BL = data_BL(:,4); %Auxillery pressure Differntial in Pa.
Eld_x_BL = data_BL(:,5); %ELD Probe x axis location in mm.
Eld_y_BL = data_BL(:,6); %ELD Probe y axis location in mm.
Voltage_BL = data_BL(:,7); % Voltage data were recorded at (in Volts).

%% Dynamic naming

% get how many different voltages we have in each file, and where they at.
%place holders


[ numvolt_VV location_VV ] = unique(Voltage_VV);
[ numvolt_BL location_BL ] = unique(Voltage_BL);

%Dynamic naming

% loop over the first file and make the names depend on the parsing
% location

set_of_names_VV = {};
set_of_names_BL = {};

for i = 1:length(numvolt_VV)
    
 set_of_names_VV{i} = [num2str(numvolt_VV(i)) '_VV' ];

end


%loop over the second file

for i = 1:length(numvolt_BL)
    
 set_of_names_BL{i} = [num2str(numvolt_BL(i)) '_BL' ];

end


%% density values

%density for VV (Velocity Voltage and Boundary layer, each will be done in
%a seperate loop.


% we will store all the matrices in a giant matrix, each row represents the density values at that
% voltage.

Density_VV = zeros(length(Voltage_VV),1);
Density_BL = zeros(length(Voltage_BL),1);


% we will run a loop to compute the densities and place them

% each input file will have its own loop


% the vleocity voltage loop
for i=1:length(numvolt_VV)
    
    %condition for the last itteration
    if i==length(numvolt_VV)
        
        for j = location_VV(i):length(Voltage_VV)
        Density_VV(j) = atm_P_VV(j)/( atm_Temp_VV(j) * RAir) ;
        end
        
    else
        
    for j = location_VV(i):location_VV(i+1)
        
        Density_VV(j) = atm_P_VV(j)/( atm_Temp_VV(j) * RAir) ;
    end
    
    end
    
end


% the boundary layer loop

for i=1:length(numvolt_BL)
    
    %condition for the last itteration
    if i==length(numvolt_BL)
        
        for j = location_BL(i):length(Voltage_BL)
        Density_BL(j) = atm_P_BL(j)/( atm_Temp_BL(j) * RAir) ;
        end
        
    else
        
    for j = location_BL(i):location_BL(i+1)
        
        Density_BL(j) = atm_P_BL(j)/( atm_Temp_BL(j) * RAir) ;
    end
    
    end
    
end





%% get the mean values: VV


% STARTING WITH VV!

% mean values for atm pressure and std.
mean_atm_P_VV_1 = mean(atm_P_VV( 1:location_VV(2)));
mean_atm_P_VV_3 = mean(atm_P_VV(location_VV(2):location_VV(3)));
mean_atm_P_VV_5 = mean(atm_P_VV(location_VV(3):location_VV(4)));
mean_atm_P_VV_7 = mean(atm_P_VV(location_VV(4):location_VV(5)));
mean_atm_P_VV_9 = mean(atm_P_VV(location_VV(5):length(atm_P_VV)));


std_atm_P_VV_1 = std(atm_P_VV( 1:location_VV(2)));
std_atm_P_VV_3 = std(atm_P_VV(location_VV(2):location_VV(3)));
std_atm_P_VV_5 = std(atm_P_VV(location_VV(3):location_VV(4)));
std_atm_P_VV_7 = std(atm_P_VV(location_VV(4):location_VV(5)));
std_atm_P_VV_9 = std(atm_P_VV(location_VV(5):length(atm_P_VV)));

% temp

mean_atm_Temp_VV_1 = mean(atm_Temp_VV( 1:location_VV(2)));
mean_atm_Temp_VV_3 = mean(atm_Temp_VV(location_VV(2):location_VV(3)));
mean_atm_Temp_VV_5 = mean(atm_Temp_VV(location_VV(3):location_VV(4)));
mean_atm_Temp_VV_7 = mean(atm_Temp_VV(location_VV(4):location_VV(5)));
mean_atm_Temp_VV_9 = mean(atm_Temp_VV(location_VV(5):length(atm_Temp_VV)));

std_atm_Temp_VV_1 = std(atm_Temp_VV( 1:location_VV(2)));
std_atm_Temp_VV_3 = std(atm_Temp_VV(location_VV(2):location_VV(3)));
std_atm_Temp_VV_5 = std(atm_Temp_VV(location_VV(3):location_VV(4)));
std_atm_Temp_VV_7 = std(atm_Temp_VV(location_VV(4):location_VV(5)));
std_atm_Temp_VV_9 = std(atm_Temp_VV(location_VV(5):length(atm_Temp_VV)));


% aux diff

mean_Aux_diff_P_VV_1 = mean(Aux_diff_P_VV( 1:location_VV(2)));
mean_Aux_diff_P_VV_3 = mean(Aux_diff_P_VV(location_VV(2):location_VV(3)));
mean_Aux_diff_P_VV_5 = mean(Aux_diff_P_VV(location_VV(3):location_VV(4)));
mean_Aux_diff_P_VV_7 = mean(Aux_diff_P_VV(location_VV(4):location_VV(5)));
mean_Aux_diff_P_VV_9 = mean(Aux_diff_P_VV(location_VV(5):length(Aux_diff_P_VV)));

std_Aux_diff_P_VV_1 = std(Aux_diff_P_VV( 1:location_VV(2)));
std_Aux_diff_P_VV_3 = std(Aux_diff_P_VV(location_VV(2):location_VV(3)));
std_Aux_diff_P_VV_5 = std(Aux_diff_P_VV(location_VV(3):location_VV(4)));
std_Aux_diff_P_VV_7 = std(Aux_diff_P_VV(location_VV(4):location_VV(5)));
std_Aux_diff_P_VV_9 = std(Aux_diff_P_VV(location_VV(5):length(Aux_diff_P_VV)));


% air diff

mean_air_diff_P_VV_1 = mean(air_diff_P_VV( 1:location_VV(2)));
mean_air_diff_P_VV_3 = mean(air_diff_P_VV(location_VV(2):location_VV(3)));
mean_air_diff_P_VV_5 = mean(air_diff_P_VV(location_VV(3):location_VV(4)));
mean_air_diff_P_VV_7 = mean(air_diff_P_VV(location_VV(4):location_VV(5)));
mean_air_diff_P_VV_9 = mean(air_diff_P_VV(location_VV(5):length(air_diff_P_VV)));

std_air_diff_P_VV_1 = std(air_diff_P_VV( 1:location_VV(2)));
std_air_diff_P_VV_3 = std(air_diff_P_VV(location_VV(2):location_VV(3)));
std_air_diff_P_VV_5 = std(air_diff_P_VV(location_VV(3):location_VV(4)));
std_air_diff_P_VV_7 = std(air_diff_P_VV(location_VV(4):location_VV(5)));
std_air_diff_P_VV_9 = std(air_diff_P_VV(location_VV(5):length(air_diff_P_VV)));



%% Getting the mean values, BL

% STARTING WITH BL!

% mean values for atm pressure and std.
mean_atm_P_BL_5 = mean(atm_P_BL(location_BL(1):length(atm_P_BL)));
std_atm_P_BL_5 = std(atm_P_BL(location_BL(1):length(atm_P_BL)));

% temp

mean_atm_Temp_BL_5 = mean(atm_Temp_BL(location_BL(1):length(atm_Temp_BL)));
std_atm_Temp_BL_5 = std(atm_Temp_BL(location_BL(1):length(atm_Temp_BL)));


% aux diff
mean_Aux_diff_P_BL_5 = mean(Aux_diff_P_BL(location_BL(1):length(Aux_diff_P_BL)));
std_Aux_diff_P_BL_5 = std(Aux_diff_P_BL(location_BL(1):length(Aux_diff_P_BL)));


% air diff

mean_air_diff_P_BL_5 = mean(air_diff_P_BL(location_BL(1):length(air_diff_P_BL)));
std_air_diff_P_BL_5 = std(air_diff_P_BL(location_BL(1):length(air_diff_P_BL)));


%% Velocity: pitot-static probe

% get the velocities for different voltages for Pitot-static and 5V for
% Boundary Later 

Velocity_VV_1_Pito = (( 2 * mean_air_diff_P_VV_1*RAir*mean_atm_Temp_VV_1) / ( mean_atm_P_VV_1))^(1/2);
Velocity_VV_3_Pito = (( 2 * mean_air_diff_P_VV_3*RAir*mean_atm_Temp_VV_3) / (( mean_atm_P_VV_3)))^(1/2) ;
Velocity_VV_5_Pito = (( 2 * mean_air_diff_P_VV_5*RAir*mean_atm_Temp_VV_5) / ( mean_atm_P_VV_5))^(1/2) ;
Velocity_VV_7_Pito = (( 2 * mean_air_diff_P_VV_7*RAir*mean_atm_Temp_VV_7) / ( mean_atm_P_VV_7))^(1/2 );
Velocity_VV_9_Pito = (( 2 * mean_air_diff_P_VV_9*RAir*mean_atm_Temp_VV_9) / ( mean_atm_P_VV_9))^(1/2 );

Velocity_BL_5_Pito = (( 2 * mean_air_diff_P_BL_5*RAir*mean_atm_Temp_BL_5) / ( mean_atm_P_BL_5))^(1/2) ;

%% Velocity: Ventori tube
AreaRatio = 1/9.5 ; 
diffP_Reading = [ 0.05 ; 0.42 ; 1.5 ; 2.9 ; 4.9 ] ;
diffP_Uncertainty = [ 0.01 ; 0.05 ; 0.05 ; 0.05 ;0.05 ];
% convert inches of water to Pascal
diffP_Reading = diffP_Reading .* 248.84;
diffP_Uncertainty = diffP_Uncertainty .* 248.84;

%get the velocity readings
Velocity_VV_1_Vento = ((( 2 * diffP_Reading(1)*RAir*mean_atm_Temp_VV_1)) / ( mean_atm_P_VV_1* ( 1- (AreaRatio)^2))) ^(1/2) ;
Velocity_VV_3_Vento = (( 2 * diffP_Reading(2)*RAir*mean_atm_Temp_VV_3) / ( mean_atm_P_VV_3* ( 1- (AreaRatio)^2))) ^(1/2) ;
Velocity_VV_5_Vento = (( 2 * diffP_Reading(3)*RAir*mean_atm_Temp_VV_5) / ( mean_atm_P_VV_5* ( 1- (AreaRatio)^2))) ^(1/2) ;
Velocity_VV_7_Vento = (( 2 * diffP_Reading(4)*RAir*mean_atm_Temp_VV_7) / ( mean_atm_P_VV_7* ( 1- (AreaRatio)^2))) ^(1/2) ;
Velocity_VV_9_Vento = (( 2 * diffP_Reading(5)*RAir*mean_atm_Temp_VV_9) / ( mean_atm_P_VV_9* ( 1- (AreaRatio)^2))) ^(1/2) ;
    
%% error calculations: venturi tube

% venturi tube: this's the resultant error equation
% the instructiosn says to use 0.25 as error for temp and ignore the
% reading uncertinity and only count the systmatic error.

error_Vento_1 = (Velocity_VV_1_Vento/2) * ( sqrt ( (SigmaatmPressure/mean_atm_P_VV_1)^2 + ( SigmaTemp/mean_atm_Temp_VV_1)^2 + (SigmaManometer/diffP_Reading(1))^2 ));
error_Vento_3 = (Velocity_VV_3_Vento/2) * ( sqrt ( (SigmaatmPressure/mean_atm_P_VV_3)^2 + ( SigmaTemp/mean_atm_Temp_VV_3)^2 + (SigmaManometer/diffP_Reading(2))^2 ));
error_Vento_5 = (Velocity_VV_5_Vento/2) * ( sqrt ( (SigmaatmPressure/mean_atm_P_VV_5)^2 + ( SigmaTemp/mean_atm_Temp_VV_5)^2 + (SigmaManometer/diffP_Reading(3))^2 ));
error_Vento_7 = (Velocity_VV_7_Vento/2) * ( sqrt ( (SigmaatmPressure/mean_atm_P_VV_7)^2 + ( SigmaTemp/mean_atm_Temp_VV_7)^2 + (SigmaManometer/diffP_Reading(4))^2 ));
error_Vento_9 = (Velocity_VV_9_Vento/2) * ( sqrt ( (SigmaatmPressure/mean_atm_P_VV_9)^2 + ( SigmaTemp/mean_atm_Temp_VV_9)^2 + (SigmaManometer/diffP_Reading(5))^2 ));

%% error calculations: pitot-static tube and Boundary Layer

% error should decrease 
% error in pitot static tube and Boundary Layer since same equation was
% used to 
error_pitot_1 = ( Velocity_VV_1_Pito/2 ) * ( sqrt ( (SigmaDiffPressure/mean_air_diff_P_VV_1)^2 + (SigmaTemp/mean_atm_Temp_VV_1)^2 + (SigmaatmPressure/mean_atm_P_VV_1)^2 )) ;
error_pitot_3 = ( Velocity_VV_3_Pito/2 ) * ( sqrt ( (SigmaDiffPressure/mean_air_diff_P_VV_3)^2 + (SigmaTemp/mean_atm_Temp_VV_3)^2 + (SigmaatmPressure/mean_atm_P_VV_3)^2 )) ;
error_pitot_5 = ( Velocity_VV_5_Pito/2 ) * ( sqrt ( (SigmaDiffPressure/mean_air_diff_P_VV_5)^2 + (SigmaTemp/mean_atm_Temp_VV_5)^2 + (SigmaatmPressure/mean_atm_P_VV_5)^2 )) ;
error_pitot_7 = ( Velocity_VV_7_Pito/2 ) * ( sqrt ( (SigmaDiffPressure/mean_air_diff_P_VV_7)^2 + (SigmaTemp/mean_atm_Temp_VV_7)^2 + (SigmaatmPressure/mean_atm_P_VV_7)^2 )) ;
error_pitot_9 = ( Velocity_VV_9_Pito/2 ) * ( sqrt ( (SigmaDiffPressure/mean_air_diff_P_VV_9)^2 + (SigmaTemp/mean_atm_Temp_VV_9)^2 + (SigmaatmPressure/mean_atm_P_VV_9)^2 )) ;

error_BL_5 = ( Velocity_BL_5_Pito/2 ) * ( sqrt ( (SigmaDiffPressure/mean_air_diff_P_VV_5)^2 + (SigmaTemp/mean_atm_Temp_VV_5)^2 + (SigmaatmPressure/mean_atm_P_VV_5)^2 )) ;

%% printout the results:

Voltage = { '1';'3';'5';'7';'9'};
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


end