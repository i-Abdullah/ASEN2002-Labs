%% info
%{ 

This script is meant to analyze the data obtained from CU-Boulder ITLL Wind
Tunnel for a Clark Y airfoil, and run aerodynamic analysis. Part of Fall 18
ASEN 2002 Labs

Done By:

1- Sam D'Souza
2- Trevor Slack
3- Foster Greer
4 - Nathan Portman
5 - Abdulla Alameri


% This codes will read .csv file, and plot coefficient of drag for profile
drag and compare it with x/c (location/cord length) and also calculate
coefficients of friciton as well as some of the 

%}



%% housekeeping

clear;
clc;
close all;

%% Read Data:

Section = 13;
Group = 15;


if Group <= 9
    filename = ['AirfoilPressure_S0' num2str(Section) '_G0' num2str(Group) '.csv' ];
else
    filename = ['AirfoilPressure_S0' num2str(Section) '_G' num2str(Group) '.csv' ];
end


%filename = [ 'AirfoilPressure_S011_G15.csv'] ;
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
SP9 = Data(:,15);
SP10 = Data(:,16);
SP11 = Data(:,17); 
SP12 = Data(:,18);
SP13 = Data(:,19);
SP14 = Data(:,20);
SP15 = Data(:,21);
SP16 = Data(:,22);

SP_All = Data(:,(7:22));

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

LocationTrail = 3.5; %where the trail is in meters.

%each two conscutive ports on the same height will be in the same loop.
%so 10 and 8 together, and 12 and 14 together.

% SP 8 AND SP 10
for i=1:1:12
y = [ mean(SP8((i*20-19):i*20)); mean(SP10((i*20-19):i*20)) ];

%get max/min values from trailing edge from upper half of the airfoil. 
y_trail_max_upper = [ max(SP8((i*20-19):i*20)); max(SP10((i*20-19):i*20)) ];
y_trail_min_upper = [ min(SP8((i*20-19):i*20)); min(SP10((i*20-19):i*20)) ];

t = [ SP8_location; SP10_location ];

%for mean values
Slope = (y(2)-y(1))/(t(2)-t(1));
Intercept = y(1) - Slope*t(1) ;

%for max and min
Slope_max = (y_trail_max_upper(2)-y_trail_max_upper(1))/(t(2)-t(1));
Intercept_max = y_trail_max_upper(1) - Slope*t(1) ;

Slope_min = (y_trail_min_upper(2)-y_trail_min_upper(1))/(t(2)-t(1));
Intercept_min = y_trail_min_upper(1) - Slope*t(1) ;

% get p values from linear extrapolation.
Ptrail_UpperPorts(i) = Slope*(LocationTrail) + Intercept;
Ptrail_UpperPorts_max(i) = Slope_max*(LocationTrail) + Intercept_max;
Ptrail_UpperPorts_min(i) = Slope_min*(LocationTrail) + Intercept_min;

%Uncertinity_Ptrail_Upper(i) = sqrt( [ LocationTrail 1] * Q * [ LocationTrail ; 1 ] );

end

for i=1:1:12
y = [ mean(SP12((i*20-19):i*20)); mean(SP14((i*20-19):i*20))];
%get max/min values from trailing edge from lower half of the airfoil. 
y_trail_max_lower = [ max(SP12((i*20-19):i*20)); max(SP14((i*20-19):i*20)) ];
y_trail_min_lower = [ min(SP12((i*20-19):i*20)); min(SP14((i*20-19):i*20)) ];

t = [ SP12_location; SP14_location ];

%calculate the new pressure at the trailing edge.

Ptrail_LowerPorts(i) = Slope*(LocationTrail) + Intercept;
%Uncertinity_Ptrail_Lower(i) = sqrt( [ LocationTrail 1] * Q * [ LocationTrail ; 1 ] );

%for mean values
Slope = (y(2)-y(1))/(t(2)-t(1));
Intercept = y(1) - Slope*t(1) ;


%for max and min
Slope_max = (y_trail_max_lower(2)-y_trail_max_lower(1))/(t(2)-t(1));
Intercept_max = y_trail_max_lower(1) - Slope*t(1) ;

Slope_min = (y_trail_min_lower(2)-y_trail_min_lower(1))/(t(2)-t(1));
Intercept_min = y_trail_min_lower(1) - Slope*t(1) ;


% get p values from linear extrapolation.
Ptrail_LowerPorts(i) = Slope*(LocationTrail) + Intercept;

Ptrail_LowerPorts_max(i) = Slope_max*(LocationTrail) + Intercept_max;
Ptrail_LowerPorts_min(i) = Slope_min*(LocationTrail) + Intercept_min;

end

% get the mean of the the two values, as well as the maxes and mins.
PTrail = mean( [ Ptrail_UpperPorts ; Ptrail_LowerPorts ] );
PTrailmax = mean( [ Ptrail_UpperPorts_max ; Ptrail_LowerPorts_max ] );
PTrailmin = mean([ Ptrail_UpperPorts_min ; Ptrail_LowerPorts_min ] );

%%

%get mean data for all ports for all angel of attacks all speeds, each
%column should have 12 rows

for i=1:1:12

PortsMeans(i,:) = mean(Data ( ((i*20-19):(i*20)),(7:22) ),1);
PortsMaxes(i,:) = max(Data ( ((i*20-19):(i*20)),(7:22) ));
PortsMins(i,:) = min(Data ( ((i*20-19):(i*20)),(7:22) ));

pitotDynamicPMean(i) = mean(pitotDynamicP ( ((i*20-19):(i*20)),1 ),1);
TatmMean(i) = mean(Patm ( ((i*20-19):(i*20)),1 ),1);
PatmMean(i) = mean(Tatm ( ((i*20-19):(i*20)),1 ),1);


end


%Loop To re-create the SP Values after getting the last port:

[ r c ] = size(PortsMeans);

% the max and mins are for error bounds.
SP_All_Updated = zeros(r,c+1);
SP_All_Updated_max = zeros(r,c+1);
SP_All_Updated_min = zeros(r,c+1);

for i = 1:c+1
    if i==10
        SP_All_Updated(:,i) = PTrail;
        SP_All_Updated_max(:,i) = PTrailmax;
        SP_All_Updated_min(:,i) = PTrailmin;

    elseif i<10 % the 10th location is where the trailing edge is @.
    SP_All_Updated(:,i) = PortsMeans(:,i);
    SP_All_Updated_max(:,i) = PortsMaxes(:,i);
    SP_All_Updated_min(:,i) = PortsMins(:,i);

    else %after the 10 is added, the index will exceed the matrix deminsions
    
     SP_All_Updated(:,i) = PortsMeans(:,i-1);
     SP_All_Updated_max(:,i) = PortsMaxes(:,i-1);
     SP_All_Updated_min(:,i) = PortsMins(:,i-1);

    end
    
end


a=1;
%% Geometry of airfoil

AirfoilGeometry = xlsread('Data/AirfoilGeometry.xlsx',1);
PortsAndConnection = xlsread('Data/AirfoilGeometry.xlsx',2);
CordLength = 3.5;

x_c = PortsAndConnection(:,2)/CordLength;

%exclude not used ports:
x_c(9)= [];
%x_c(11)= [];
x_c(13)= [];
x_c(15)= [];


%Calculate ports mean values for pressure @ airfoil
for i=1:1:12
    
 %each column represent 1 air speed per angle of attack
 
 % thus each row is one port, we have 16 ports, thus 16 raws
 % 4 angel of attacks *  3 speeds, thus 12 columns.
 
 
PortsMeanValues(:,i) = mean(Data ( ((i*20-19):(i*20)),(7:22) ),1);

% Angel of attack related to each set of columns

% this will b 2 rows * N columns, the first row is angel of attack, second
% is air speed. 

Alpha_V(:,i) = [ mean(Data ( ((i*20-19):(i*20)),(23) ),1) ; mean(Data ( ((i*20-19):(i*20)),(4) ),1) ];


% get qinfinity and Pinf, descreption below.

Pinf(:,i) = mean(Data ( ((i*20-19):(i*20)),(6) ),1);
qinf(:,i) = mean(Data ( ((i*20-19):(i*20)),(5) ),1);


end

a = 1;
%Now we can calculate Cp (Coefficient of pressure)

%{

Cp = (P - Pinf)/(qinf);
P = PortsMeanValues;
Pinf = Air pressure transducer 6
Qinfity = Dynamic pressure from pitot tube (@ column 5)

%}
for i=1:1:12
    
  Cp(:,i) = (SP_All_Updated(i,:)')./(qinf(:,i));
  Cp_max(:,i) = (SP_All_Updated_max(i,:)')./(qinf(:,i));
  Cp_min(:,i) = (SP_All_Updated_min(i,:)')./(qinf(:,i));


end


%% Lift and Drag


% get ports locations
PortsExlusion = PortsAndConnection(:,:);


%exclude not used pots
PortsExlusion(9,:) = [];
PortsExlusion(13,:) = [];
PortsExlusion(15,:) = [];

dx=zeros(length(PortsExlusion),1);
dy=zeros(length(PortsExlusion),1);

for i = 1:length(PortsExlusion)-1 %calculate the delta x and delta y values 
    dx(i) = PortsExlusion(i+1,2)-PortsExlusion(i,2);
    dy(i) = PortsExlusion(i+1,3)-PortsExlusion(i,3);
end



[ r c ] = size(Cp);
for j=1:c
    for i=1:r-1
        
    Cn(i,j) = 0.5*(Cp(i,j) + Cp(i+1,j))* (dx(i)/CordLength) ;
    Ca(i,j) = 0.5*(Cp(i,j) + Cp(i+1,j))* (dy(i)/CordLength) ;
    Cn_max(i,j) = 0.5*(Cp_max(i,j) + Cp_max(i+1,j))* (dx(i)/CordLength) ;
    Cn_min(i,j) = 0.5*(Cp_min(i,j) + Cp_min(i+1,j))* (dx(i)/CordLength) ;
    Ca_max(i,j) = 0.5*(Cp_max(i,j) + Cp_max(i+1,j))* (dy(i)/CordLength) ;
    Ca_min(i,j) = 0.5*(Cp_min(i,j) + Cp_min(i+1,j))* (dy(i)/CordLength) ;

    end
end

%sum all values
Cn = -sum(Cn,1);
Ca = sum(Ca,1);

Cn_max = -sum(Cn_max,1);
Ca_max = sum(Ca_max,1);

Cn_min = -sum(Cn_min,1);
Ca_min = sum(Ca_min,1);

%coefficient of lift and drag
Cl = (Cn .* cosd(Alpha_V(1,:))) - (Ca .* sind(Alpha_V(1,:))) ;
Cd = (Ca .* cosd(Alpha_V(1,:))) + (Cn .* sind(Alpha_V(1,:))) ;

Cl_max = (Cn_max .* cosd(Alpha_V(1,:))) - (Ca_max .* sind(Alpha_V(1,:))) ;
Cd_max = (Ca_max .* cosd(Alpha_V(1,:))) + (Cn_max .* sind(Alpha_V(1,:))) ;


Cl_min = (Cn_min .* cosd(Alpha_V(1,:))) - (Ca_min .* sind(Alpha_V(1,:))) ;
Cd_min = (Ca_min .* cosd(Alpha_V(1,:))) + (Cn_min .* sind(Alpha_V(1,:))) ;

% now, we must note that the mins and maxe matrixes has the maximum values
% and minumums for all ports, and hence, for instance, @ port 1, there's a
% max and a min, the max is stored in one matrix, the min is stored in
% anther, thus 
%% plotting cp with respect to x_c


%plot the first set
for i=1:3
figure(1)
subplot(3,1,i)
plot(x_c(:,1),Cp(:,i),'*')
hold on
plot(x_c(:,1),zeros(1,length(Cp(:,i))),'k')
patch([x_c(:,1), fliplr(x_c(:,1))], [Cp(:,i) fliplr(Cp(:,i))], [0.9290, 0.6940, 0.1250], 'FaceAlpha',0.8)
alpha(0.1); %change transperancy of  the filling

hold on
%plot maximas and minmas 

%maxes and mins
errorbar(x_c(:,1),Cp(:,i),Cp_max(:,i),Cp_min(:,i))
hold on


grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
legend({'Cp Values','Linear connection','Area under the curve','Error bars'},'Orientation','horizontal','Location','SouthEast');
xlabel( ' x/c ')
ylabel ('Cp')

set(gca, 'YDir','reverse')


hold off;



end

% the second portion

for i=4:6
figure(2)
subplot(3,1,i-3)
plot(x_c(:,1),Cp(:,i),'*')
hold on
plot(x_c(:,1),zeros(1,length(Cp(:,i))),'k')
patch([x_c(:,1), fliplr(x_c(:,1))], [Cp(:,i) fliplr(Cp(:,i))], [0.9290, 0.6940, 0.1250], 'FaceAlpha',0.8)
alpha(0.1); %change transperancy of  the filling

hold on
errorbar(x_c(:,1),Cp(:,i),Cp_max(:,i),Cp_min(:,i))
hold on

grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
xlabel( ' x/c ')
legend({'Cp Values','Linear connection','Area under the curve','Error bars'},'Orientation','horizontal','Location','SouthEast');
ylabel ('Cp')

set(gca, 'YDir','reverse')

hold off;

end


%the third portion

for i=7:9
figure(3)
subplot(3,1,i-6)
plot(x_c(:,1),Cp(:,i),'*')
hold on
plot(x_c(:,1),zeros(1,length(Cp(:,i))),'k')
patch([x_c(:,1), fliplr(x_c(:,1))], [Cp(:,i) fliplr(Cp(:,i))], [0.9290, 0.6940, 0.1250], 'FaceAlpha',0.8)
alpha(0.1); %change transperancy of  the filling

hold on
errorbar(x_c(:,1),Cp(:,i),Cp_max(:,i),Cp_min(:,i))
hold on
grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
xlabel( ' x/c ')
legend({'Cp Values','Linear connection','Area under the curve','Error bars'},'Orientation','horizontal','Location','SouthEast');
ylabel ('Cp')

set(gca, 'YDir','reverse')

hold off;

end

%forth portion

for i=10:12
figure(4)
subplot(3,1,i-9)
plot(x_c(:,1),Cp(:,i),'*')
hold on
plot(x_c(:,1),zeros(1,length(Cp(:,i))),'k')
patch([x_c(:,1), fliplr(x_c(:,1))], [Cp(:,i) fliplr(Cp(:,i))], [0.9290, 0.6940, 0.1250], 'FaceAlpha',0.8)
alpha(0.1); %change transperancy of  the filling

hold on
errorbar(x_c(:,1),Cp(:,i),Cp_max(:,i),Cp_min(:,i))
hold on

grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
xlabel( ' x/c ')
legend({'Cp Values','Linear connection','Area under the curve','Error bars'},'Orientation','horizontal','Location','SouthEast');
ylabel ('Cp')

set(gca, 'YDir','reverse')

hold off;

end

%% plot alpha vs cl/cd (same like L/D;

figure(5);

scatter(Alpha_V(1,:),(Cl./Cd))
grid minor
title(['\alpha Vs L/D']);


%% error analysis: Cp: Not needed.

% uncertainties values are the same like the previous lab, check lab report
% to see how those are done.


% SigmaTemp = 0.25 ; % in k
% SigmaatmPressure = (250-20)*10^3*(1.5/100); %from lab document
% SigmaDiffPressure = 6894.76 * (1/100); %from lab document
% SgimaScav = 0.20;
% % get error fro v_inifity
% Rfluid = 287;
% 
% % Use pitot function
% [ Velocinf Errorinf ] = Pitot (PatmMean, TatmMean, pitotDynamicPMean, SigmaatmPressure, SigmaTemp, SigmaDiffPressure,Rfluid);
% 
% 
% 
% % now get the error for Cp:
% 
% 
% syms Pscav Patmo Tatmo V_freeStream
% 
% Cp_error_eqn(Pscav,Patmo,Tatmo,V_freeStream) = Pscav / ((0.5)*(Patmo/(Rfluid*Tatmo))*(V_freeStream));
% 
% % get partial derivatives
% 
% Cp_Partial_Pscav = diff(Cp_error_eqn,Pscav);
% Cp_Partial_Patmo = diff(Cp_error_eqn,Patmo);
% Cp_Partial_Tatmo = diff(Cp_error_eqn,Tatmo);
% Cp_Partial_V_freeStream = diff(Cp_error_eqn,V_freeStream);
% 
% 
% % if there's constant uncertinity values, create a matrix with their
% % size
%     
% 
% % if length(SigmaDiffPressure) == 1 && length(SigmaTemp)==1 && length(SigmaatmPressure) == 1
% %     SigmaDiffPressure = ones(1,length(P_atm)) .* SigmaDiffPressure;
% %     SigmaTemp = ones(1,length(P_atm)) .* SigmaDiffPressure;
% %     SigmaatmPressure = ones(1,length(P_atm)) .* SigmaDiffPressure;
% % 
% % else
% %     
% %     SigmaDiffPressure = sigma_Air_P_Diff;
% %     SigmaTemp = SigmaTemp;
% %     SigmaatmPressure = SigmaatmPressure;
% %     
% % end
% 
% 
% % loop to get all errors: one 
% % this will be loop for each speed/angel of attack once.
% 
% [ r c ] = size(Cp);
% for j=1:c
%     for i=1:r-1
%         
%  value_before_root(i,j) = ( ((Cp_Partial_Pscav( SP_All_Updated(j,i) ,PatmMean(j),TatmMean(j),Velocinf(j)))*SgimaScav)^2 ...
%  + ((Cp_Partial_Patmo(SP_All_Updated(j,i),PatmMean(j),TatmMean(j),Velocinf(j)))*SigmaatmPressure)^2 ...
%  +((Cp_Partial_Tatmo(SP_All_Updated(j,i),PatmMean(j),TatmMean(j),Velocinf(j)))*SigmaTemp)^2 ...
%  +((Cp_Partial_V_freeStream(SP_All_Updated(j,i),PatmMean(j),TatmMean(j),Velocinf(j)))*Errorinf(j))^2 ) ;
%  
%  Erros_Cp(i,j) = double(sqrt(value_before_root(i,j)));
% 
%     end
% end
% 
% % now we have to do the last port error, which is combinations of two
% % edges
% 
% % the Error_Cp has 16*12, each column is one test, each row is one of the
% % ports in order.
% 
% % the trailing edge has 
% 
