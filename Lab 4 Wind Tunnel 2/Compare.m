function [Cd,Cl,alpha]=Compare(Section,Group)
% This function is a the exact same script WindTunnel.m but turned into
% a function that runs by the group and section number, to analyze a .csv
% file that is obtained from CU Wind Tunnel experiments. 
%
% ------------------------ (INPUTS) --------------
%
%       1- Section Number: Between 11 and 14) 
%
%       2- Group Number: Odd numbers, starting from 1, ending at
%       15 for each section
%
% ------------------------ (OUTPUTS) --------------
%
%       1- Cd: coefficient of drag, an array that has the size of 1x(n*k)
%              where n is the number of angel of attacks and k is the
%              number of speeds tested, thus in our case here it is 1x12,
%              where say for instance we have three angel of attacks and 4
%              different speeds, the first three elments are for the first
%              speed, the second three are for the third, and so on.
%
%
%       2- Cl: coefficient of lift, an array that has the size of 1x(n*k)
%              where n is the number of angel of attacks and k is the
%              number of speeds tested, thus in our case here it is 1x12,
%              where say for instance we have three angel of attacks and 4
%              different speeds, the first three elments are for the first
%              speed, the second three are for the third, and so on.
%
%       3- alpha: angel of attack, it will also have the same sizing, and
%       thus in this case the first three elements will represnt the first
%       the same angel of attack because the only thing that changes is
%       speed, and so on.
%
%
% ------------------- ( Done By )------------------
%           1- Sam D'Souza
%           2- Trevor Slack
%           3- Foster Greer
%           4 - Nathan Portman
%           5 - Abdulla Alameri


%% Read Data:


if Group <= 9
    filename = ['AirfoilPressure_S0' num2str(Section) '_G0' num2str(Group) '.csv' ];
    if isfile(['Data/' filename ]) ~= 1
        Cd=0;
        Cl=0;
        alpha=0;
        return
    end
else
    filename = ['AirfoilPressure_S0' num2str(Section) '_G' num2str(Group) '.csv' ];
    if isfile(['Data/' filename ]) ~= 1
        Cd=0;
        Cl=0;
        alpha=0;
        return
    end
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
t = [ SP8_location; SP10_location ];

Slope = (y(2)-y(1))/(t(2)-t(1));
Intercept = y(1) - Slope*t(1) ;

Ptrail_UpperPorts(i) = Slope*(LocationTrail) + Intercept;
%Uncertinity_Ptrail_Upper(i) = sqrt( [ LocationTrail 1] * Q * [ LocationTrail ; 1 ] );

end

for i=1:1:12
y = [ mean(SP12((i*20-19):i*20)); mean(SP14((i*20-19):i*20)) ];
t = [ SP12_location; SP14_location ];

%calculate the new pressure at the trailing edge.

Ptrail_LowerPorts(i) = Slope*(LocationTrail) + Intercept;
%Uncertinity_Ptrail_Lower(i) = sqrt( [ LocationTrail 1] * Q * [ LocationTrail ; 1 ] );

Slope = (y(2)-y(1))/(t(2)-t(1));
Intercept = y(1) - Slope*t(1) ;

Ptrail_LowerPorts(i) = Slope*(LocationTrail) + Intercept;

end

PTrail = mean( [ Ptrail_LowerPorts ; Ptrail_LowerPorts ] );


%get mean data for all ports for all angel of attacks all speeds, each
%column should have 12 rows

for i=1:1:12

PortsMeans(i,:) = mean(Data ( ((i*20-19):(i*20)),(7:22) ),1);

end


%Loop To re-create the SP Values after getting the last port:

[ r c ] = size(PortsMeans);

SP_All_Updated = zeros(r,c+1);
for i = 1:c+1
    if i==10
        SP_All_Updated(:,i) = PTrail;
    elseif i<10
    SP_All_Updated(:,i) = PortsMeans(:,i);
    
    else %after the 10 is added, the index will exceed the matrix deminsions
    
     SP_All_Updated(:,i) = PortsMeans(:,i-1);

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
    
  Cp(:,i) = (SP_All_Updated(i,:)' - Pinf(:,i))./(qinf(:,i));

end


%% Lift and Drag

PortsExlusion = PortsAndConnection(:,:);


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
    
    end
end

%sum all values
Cn = -sum(Cn,1);
Ca = sum(Ca,1);

%coefficient of lift and drag
Cl = (Cn .* cosd(Alpha_V(1,:))) - (Ca .* sind(Alpha_V(1,:))) ;
Cd = (Ca .* cosd(Alpha_V(1,:))) + (Cn .* sind(Alpha_V(1,:))) ;


%% store Alpha values.

alpha=Alpha_V(1,:);
end