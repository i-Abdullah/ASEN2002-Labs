function [ Results, BLThickness, ErrorInBLThickness ] = Lab3(SectionNum,GroupNum)
%% info
%{

Thie function will extract and analyze data obtained from a wind tunnel lab,
part of ASEN 2002: Lab 3, CU Boulder, Fall 18. It is designed to read csv
files from class of 2021 CU Boulder. The files are stored in a folder
called Data and named the same with different group/section numbers, and
thus all you need to do here is to give section number and group number and
the function will do the rest. Some adjustments needs to be made manually
for the errors, manometer readings if used, and what goes inside each
function aroun line 163,164.


Done by

1- Jack Soltys
2- Abdulla Al Ameri
3- Greer Foster
4- Caelan Maitland


Input:
1- Section number
2- Group Number

Output:

Results : Table of voltage vs velocities and errors.
BLThickness : Boundary Layer thickness in units of y-probe location in file input
ErrorInBLThickness : Error in Boundary Layer thickness.

%}


%% define constants/ hard code

% those are constants, change them here.


RAir = 287.0 ; % pa / m^3 k
SigmaTemp = 0.25 ; % in k
SigmaatmPressure = (250-20)*10^3*(1.5/100); %from lab document
SigmaDiffPressure = 6894.76 * (1/100); %from lab document
SigmaManometer = 0.1*248 ; % in inch
AreaRatio = 1/9.5 ;


%manometer readings

%248 to convert to Pascal from inches of water.

ManoReadings = [ 0.05 ; 0.42 ; 1.5 ; 2.9 ; 4.9 ] .* 248;
ManoReadings = [ 0.01 ; 0.25 ; 1.15 ; 2.54 ; 4.4 ] .* 248;

ManoUncert = [ 0.01 ; 0.05 ; 0.05 ; 0.05 ; 0.05 ] .* 248 ;

% input files

if GroupNum <= 9
    inputfileVV = ['VelocityVoltage_S0' num2str(SectionNum) '_G0' num2str(GroupNum) '.csv' ];
    inputfileBL = ['BoundaryLayer_S0' num2str(SectionNum) '_G0' num2str(GroupNum) '.csv' ];
else
    inputfileVV = ['VelocityVoltage_S0' num2str(SectionNum) '_G' num2str(GroupNum) '.csv' ];
    inputfileBL = ['BoundaryLayer_S0' num2str(SectionNum) '_G' num2str(GroupNum) '.csv' ];
end


% title of the BL plot.
Title = [ 'Velocity vs Y-probe Location for Group: ' num2str(GroupNum) ', Section: ' num2str(SectionNum) ] ;


%% read the files/ call the reading fucntion : VV

[ VV_Files BL_Files ] = SubsonicTunnel(inputfileVV,inputfileBL,RAir);

%% Seperate 
%extract values:
Data_VV = VV_Files(:,2);;
VV_Patm = Data_VV{1,1}{:,1} ;
VV_Air_P_diff = Data_VV{2,1}{:,1};
VV_Aux_P_diff = Data_VV{3,1}{:,1};
VV_Temp = Data_VV{4,1}{:,1};
VV_Density = Data_VV{5,1}{:,1};
VV_Voltages = Data_VV{6,1}{:,1};


% put all the mean values together, all the std values together.

% we can do a trick here, since the number of voltages is the same, getting
% then all the pressure, And temp values will be the same size, thus we can
% use one loop.

[ rows columns ] = size(VV_Patm);

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

%% read the files/ call the reading fucntion : BL


%extract values:
Data_BL = BL_Files(:,2);;
BL_Ylocation = Data_BL{1,1}{:,1} ;
BL_Ylocation_std = Data_BL{2,1}{:,1};
BL_Velocity = Data_BL{3,1}{:,1};


% since all of them the same length we can use the same length

[ rows columns ] = size(BL_Ylocation);

Ylocation_BL_values = zeros(rows,1);
Ylocation_std_BL_values = zeros(rows,1);
Velocity_BL_values = zeros(rows,1);


for i=1:rows
        
Ylocation_BL_values(i) = cell2mat(BL_Ylocation(i));
Ylocation_std_BL_values(i) = cell2mat(BL_Ylocation_std(i));
Velocity_BL_values(i) = cell2mat(BL_Velocity(i));

end

%% asume constant uncertinity:

% this section can be omitted and std values should be replaced for uncertainity.


sigma_T_atm = ones(1,length(atmTemp_MeanValues_VV)) * SigmaTemp ;
sigma_Air_P_Diff = ones(1,length(Air_P_diff_MeanValues_VV)) * SigmaDiffPressure ;
sigma_P_atm = ones(1,length(Patm_MeanValues_VV)) * SigmaatmPressure ;
sigma_Manometer = ones(1,length(Patm_MeanValues_VV)) * SigmaManometer;

%% calculate Velocity


[ Velc_Venturi Error_Venturi ] = Venturi (Patm_MeanValues_VV, atmTemp_MeanValues_VV, ManoReadings, sigma_P_atm, sigma_T_atm, sigma_Manometer,RAir,AreaRatio);
[ Velc_Pitot Error_Pitot ] = Pitot (Patm_MeanValues_VV, atmTemp_MeanValues_VV, Air_P_diff_MeanValues_VV, sigma_P_atm, sigma_T_atm,sigma_Air_P_Diff,RAir);


%% Boundary Layer thickness

% make a fit, plug the 95% free stream velocity value, and get the BL
% Thickness


VelocAtBL = Velocity_BL_values(12)*0.95;
[ Function ErrorStruct ] = createFit(Velocity_BL_values, Ylocation_BL_values,Title);
hold on
plot(VelocAtBL, feval(Function,VelocAtBL), '*')
hold on
refline(0,feval(Function,VelocAtBL));
legend( 'Data points', 'Excluded Local V_\infty', 'Best Fit Line','95% Of Local V_\infty ', 'Corresponding BL thickness','Location', 'NorthWest' );
hold off

BLThickness = feval(Function,VelocAtBL);
RisdualSum = getfield(ErrorStruct,'sse'); % Risdual sum
ErrorInBLThickness = sqrt(1/9) * RisdualSum;  % sqrt of 1/N-2



figure(2)
hold on
errorbar([str2num(cell2mat(VV_Files{6,2}{:,1}))], Velc_Venturi, Error_Venturi,'-*')
errorbar([str2num(cell2mat(VV_Files{6,2}{:,1}))], Velc_Pitot, Error_Pitot,'-*')
xlim([0.5 10.5])
legend('Venturi Tube', 'Pitot-Static Probe')
title('Airspeed Measurements for Venturi Tube and Pitot-Static Probe')
xlabel('Voltage (Volts)')
ylabel('Velocity m/s')
grid minor
 
 % Use a linear fit to create a model for airspeed at a given voltage
 
x = [Velc_Venturi'];
b = [1,1;1,3;1,5;1,7;1,9];
m = b\x; % slope and y-intercept for linear model
voltagevec = 1:0.1:9;
plot(voltagevec,m(2)*voltagevec+m(1))
legend('Venturi Tube', 'Pitot-Static Probe', 'Linear Fit Model','Location','SouthEast')
hold off


%% printout the results:

Voltage = VV_Files{6,2}{:,1};
Error_Vento = { Error_Venturi(1) ; Error_Venturi(2) ; Error_Venturi(3) ; Error_Venturi(4) ; Error_Venturi(5) };
Veloc_Pitot = { Velc_Pitot(1) ; Velc_Pitot(2) ; Velc_Pitot(3) ; Velc_Pitot(4) ; Velc_Pitot(5) };
Error_Pitot = { Error_Pitot(1) ; Error_Pitot(2) ; Error_Pitot(3) ; Error_Pitot(4) ; Error_Pitot(5) };
Veloc_Venturi = { Velc_Venturi(1) ; Velc_Venturi(2) ; Velc_Venturi(3) ; Velc_Venturi(4) ; Velc_Venturi(5) };

Results = table(Voltage,Veloc_Pitot,Error_Pitot,Veloc_Venturi,Error_Vento);

%% Boundary Layer:

fprintf('Boundary Layer  thickness is: %f ', BLThickness);
fprintf( char(177) );
fprintf( ' %f ', ErrorInBLThickness ) ;
fprintf( 'mm \n' );
%}
end