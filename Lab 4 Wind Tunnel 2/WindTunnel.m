%% housekeeping

clear;
clc;
close all;

%% Read Data:

Section = 11;
Group = 15;

filename = [ 'AirfoilPressure_S011_G15.csv'] ;
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



%% Geometry of airfoil

AirfoilGeometry = xlsread('Data/AirfoilGeometry.xlsx',1);
PortsAndConnection = xlsread('Data/AirfoilGeometry.xlsx',2);
CordLength = 3.5;

x_c = PortsAndConnection(:,2)/CordLength;

%exclude not used ports:
x_c(9)= [];
x_c(11)= [];
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
    
  Cp(:,i) = (PortsMeanValues(:,i) - Pinf(:,i))./(qinf(:,i));

end

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

grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
xlabel( ' x/c ')
ylabel ('Cp')

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

grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
xlabel( ' x/c ')
ylabel ('Cp')

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

grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
xlabel( ' x/c ')
ylabel ('Cp')

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

grid minor
title(['\alpha = ' num2str(Alpha_V(1,i)) ', V = ' num2str(Alpha_V(2,i))]);
xlabel( ' x/c ')
ylabel ('Cp')

hold off;

end

a = 1;
