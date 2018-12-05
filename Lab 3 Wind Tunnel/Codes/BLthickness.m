%% info

% this script will go and get all boundary layer thicknesses for the lab,
% if not sure what's going on or how this's going, read the comments in
% Lab3.m. Furthermore, the code will plot and compare the thicknesses.

%% housekeeping

clear;
clc;
close all;

%% get all the bondary layer thicknesses

% ther's 4 sections, 11,12,13,14, and the group numbers are odd numbers,
% 1,3,5,7.........15

groups = 1:2:15;
BL_thickness = zeros(length(groups),4);
BL_error = zeros(length(groups),4);
Vinf = zeros(length(groups),4);

for i = 11:14
    
    
    for j= 1:length(groups);
        
        %the following if statments are done to exclude some groups @
        %section i group j, because their data has errors or missing
        %important components. 
        
        
        if i==12 && groups(j)== 7;
           
            
        elseif i==13 && groups(j)==11;
            
        elseif i==14 && groups(j)==7;
      
        elseif i==12 && groups(j) == 13;
            
        elseif i==14 && groups(j) == 5;
                
        else
            
       [ ignored, Thickness, Error, Veinf ]  = Lab3Modified(i,groups(j));
        BL_thickness(j,i-10) = Thickness;
        BL_error(j,i-10) = Error;
        Vinf(j,i-10) = Veinf;


        end            
  
    
    end
    
end

% there was total 11 ports, more than one group connected to each port.

% weighted avg

%{
%% port 1

%weighted avg
 WeightedAvg_p1 = 0;
 errorWeighted_p1 = 0;
 
for i=1:2
   Weight = ( 1/(BL_error(i,1)^2) );
   WeightedAvg_p1 = WeightedAvg_p1 + ( (Weight*BL_thickness(i,1))) ;
   errorWeighted_p1 = errorWeighted_p1 + Weight;
end

WeightedAvg_p1 = WeightedAvg_p1 /errorWeighted_p1;
errorWeighted_p1 = 1 / sqrt((errorWeighted_p1));


%% port 2


%weighted avg
 WeightedAvg_p2 = 0;
 errorWeighted_p2 = 0;
 
for i=1:2
   Weight = ( 1/(BL_error(i,2)^2) );
   WeightedAvg_p2 = WeightedAvg_p2 + ( (Weight*BL_thickness(i,2))) ;
   errorWeighted_p2 = errorWeighted_p2 + Weight;
end

WeightedAvg_p2 = WeightedAvg_p2 /errorWeighted_p2;
errorWeighted_p2 = 2 / sqrt((errorWeighted_p2));


%% port 3

 WeightedAvg_p3 = 0;
 errorWeighted_p3 = 0;
 
for i=1:2
   Weight = ( 1/(BL_error(i,3)^2) );
   WeightedAvg_p3 = WeightedAvg_p3 + ( (Weight*BL_thickness(i,3))) ;
   errorWeighted_p3 = errorWeighted_p3 + Weight;
end

WeightedAvg_p3 = WeightedAvg_p3 /errorWeighted_p3;
errorWeighted_p3 = 2 / sqrt((errorWeighted_p3));


%% port 4

 WeightedAvg_p4 = 0;
 errorWeighted_p4 = 0;
 
for i=3:4
   Weight = ( 1/(BL_error(i,1)^2) );
   WeightedAvg_p4 = WeightedAvg_p4 + ( (Weight*BL_thickness(i,1))) ;
   errorWeighted_p4 = errorWeighted_p4 + Weight;
end

for i=7:8
   Weight = ( 1/(BL_error(i,4)^2) );
   WeightedAvg_p4 = WeightedAvg_p4 + ( (Weight*BL_thickness(i,4))) ;
   errorWeighted_p4 = errorWeighted_p4 + Weight;
end


WeightedAvg_p4 = WeightedAvg_p4 /errorWeighted_p4;
errorWeighted_p4 = 2 / sqrt((errorWeighted_p4));

%% port 5

 WeightedAvg_p5 = 0;
 errorWeighted_p5 = 0;
 
 %group 7 section 12 exluded, here we only have group 5 section 12 data.
for i=3
   Weight = ( 1/(BL_error(i,2)^2) );
   WeightedAvg_p5 = WeightedAvg_p5 + ( (Weight*BL_thickness(i,2))) ;
   errorWeighted_p5 = errorWeighted_p5 + Weight;
end



WeightedAvg_p5 = WeightedAvg_p5 /errorWeighted_p5;
errorWeighted_p5 = 2 / sqrt((errorWeighted_p5));

%% port 6

 WeightedAvg_p6 = 0;
 errorWeighted_p6 = 0;
 
for i=3:4
   Weight = ( 1/(BL_error(i,3)^2) );
   WeightedAvg_p6 = WeightedAvg_p6 + ( (Weight*BL_thickness(i,3))) ;
   errorWeighted_p6 = errorWeighted_p6 + Weight;
end

for i=1:2
   Weight = ( 1/(BL_error(i,4)^2) );
   WeightedAvg_p6 = WeightedAvg_p6 + ( (Weight*BL_thickness(i,4))) ;
   errorWeighted_p6 = errorWeighted_p6 + Weight;
end


WeightedAvg_p6 = WeightedAvg_p6 /errorWeighted_p6;
errorWeighted_p6 = 2 / sqrt((errorWeighted_p6));

%% port 7

 WeightedAvg_p7 = 0;
 errorWeighted_p7 = 0;
 
for i=5:6
   Weight = ( 1/(BL_error(i,1)^2) );
   WeightedAvg_p7 = WeightedAvg_p7 + ( (Weight*BL_thickness(i,1))) ;
   errorWeighted_p7 = errorWeighted_p7 + Weight;
end

for i=7:8
   Weight = ( 1/(BL_error(i,3)^2) );
   WeightedAvg_p7 = WeightedAvg_p7 + ( (Weight*BL_thickness(i,3))) ;
   errorWeighted_p7 = errorWeighted_p7 + Weight;
end


WeightedAvg_p7 = WeightedAvg_p7 /errorWeighted_p7;
errorWeighted_p7 = 2 / sqrt((errorWeighted_p7));


%% port 8

 WeightedAvg_p8 = 0;
 errorWeighted_p8 = 0;
 
for i=5:6
   Weight = ( 1/(BL_error(i,2)^2) );
   WeightedAvg_p8 = WeightedAvg_p8 + ( (Weight*BL_thickness(i,2))) ;
   errorWeighted_p8 = errorWeighted_p8 + Weight;
end

%{
%this loop is omitted because both groups group 5 and 7 from section 14 did
%not provide good data.
for i=3:4
   Weight = ( 1/(BL_error(i,4)^2) );
   WeightedAvg_p8 = WeightedAvg_p8 + ( (Weight*BL_thickness(i,4))) ;
   errorWeighted_p8 = errorWeighted_p8 + Weight;
end
%}

WeightedAvg_p8 = WeightedAvg_p8 /errorWeighted_p8;
errorWeighted_p8 = 2 / sqrt((errorWeighted_p8));


%% port 9:

 WeightedAvg_p9 = 0;
 errorWeighted_p9 = 0;
 
 %group 11 is excluded.
for i=5
   Weight = ( 1/(BL_error(i,3)^2) );
   WeightedAvg_p9 = WeightedAvg_p9 + ( (Weight*BL_thickness(i,3))) ;
   errorWeighted_p9 = errorWeighted_p9 + Weight;
end

WeightedAvg_p9 = WeightedAvg_p9 /errorWeighted_p9;
errorWeighted_p9 = 2 / sqrt((errorWeighted_p9));

%% port 10


 WeightedAvg_p10 = 0;
 errorWeighted_p10 = 0;
 
for i=7:8
   Weight = ( 1/(BL_error(i,1)^2) );
   WeightedAvg_p10 = WeightedAvg_p10 + ( (Weight*BL_thickness(i,1))) ;
   errorWeighted_p10 = errorWeighted_p10 + Weight;
end

for i=5:6
   Weight = ( 1/(BL_error(i,4)^2) );
   WeightedAvg_p10 = WeightedAvg_p10 + ( (Weight*BL_thickness(i,4))) ;
   errorWeighted_p10 = errorWeighted_p10 + Weight;
end


WeightedAvg_p10 = WeightedAvg_p10 /errorWeighted_p10;
errorWeighted_p10 = 2 / sqrt((errorWeighted_p10));



%% port 11


 WeightedAvg_p11 = 0;
 errorWeighted_p11 = 0;
 
 %section 12 group 13 is excluded.
for i=8
   Weight = ( 1/(BL_error(i,2)^2) );
   WeightedAvg_p11 = WeightedAvg_p11 + ( (Weight*BL_thickness(i,2))) ;
   errorWeighted_p11 = errorWeighted_p11 + Weight;
end



WeightedAvg_p11 = WeightedAvg_p11 /errorWeighted_p11;
errorWeighted_p11 = 2 / sqrt((errorWeighted_p11));

%% all the ports:

PortLocation = [ 9.05 10.03 11.01 11.99 12.97 13.95 14.93 15.91 16.89 17.87 18.85 ] %in inches
BLValues = [ WeightedAvg_p1 WeightedAvg_p2 WeightedAvg_p3 WeightedAvg_p4 WeightedAvg_p5 WeightedAvg_p6 WeightedAvg_p7 WeightedAvg_p8 WeightedAvg_p9 WeightedAvg_p10 WeightedAvg_p11 ]; %in mm
BLUncertinity = [ errorWeighted_p1 errorWeighted_p2 errorWeighted_p3 errorWeighted_p4 errorWeighted_p5 errorWeighted_p6 errorWeighted_p7 errorWeighted_p8 errorWeighted_p9 errorWeighted_p10 errorWeighted_p11 ];

%}


%% port 1
GroupsPerPort = 2;
%weighted avg
 Port1 = 0;
 ErrorPort1 = 0;
 Vinf1 = 0;
 
 Port1(1) = BL_thickness(1,1);
 Port1(2) = BL_thickness(2,1);
 ErrorPort1(1) = BL_error(1,1);
 ErrorPort1(2) = BL_error(2,1);
 
 Vinf1(1) = Vinf(1,1);
 Vinf1(2) = Vinf(2,1);

 YMatrix_Port1 = ones(GroupsPerPort,1)*1;

%% port 2

GroupsPerPort = 2;

%weighted avg
 Port2 = 0;
 ErrorPort2 = 0;
 
 Port2(1) = BL_thickness(1,2);
 Port2(2) = BL_thickness(2,2);
 ErrorPort2(1) = BL_error(1,2);
 ErrorPort2(2) = BL_error(2,2);
 Vinf2 = Vinf(1,2);
 Vinf2 = Vinf(2,2);

YMatrix_Port2 = ones(GroupsPerPort,1)*2;

%% port 3

GroupsPerPort = 2;

 Port3 = 0;
 ErrorPort3 = 0;
 
 Port3(1) = BL_thickness(1,3);
 Port3(2) = BL_thickness(2,3);
 ErrorPort3(1) = BL_error(1,3);
 ErrorPort3(2) = BL_error(2,3);

 Vinf3 = Vinf(1,3);
 Vinf3 = Vinf(2,3);

YMatrix_Port3 = ones(GroupsPerPort,1)*3;

%% port 4

GroupsPerPort = 4;
 Port4 = 0;
 ErrorPort4 = 0;
 
 Port4(1) = BL_thickness(3,1);
 Port4(2) = BL_thickness(4,1);
 Port4(3) = BL_thickness(7,4);
 Port4(4) = BL_thickness(8,4);
 
 ErrorPort4(1) = BL_error(3,1);
 ErrorPort4(2) = BL_error(4,1);
 ErrorPort4(3) = BL_error(7,4);
 ErrorPort4(4) = BL_error(8,4);
 
 Vinf4 = Vinf(3,1);
 Vinf4 = Vinf(4,1);
 Vinf4 = Vinf(7,4);
 Vinf4 = Vinf(8,4);

 
YMatrix_Port4 = ones(GroupsPerPort,1)*4;
%% port 5

GroupsPerPort = 1;
 Port5 = 0;
 ErrorPort5 = 0;
 
 Port5(1) = BL_thickness(3,2);
 ErrorPort5(1) = BL_error(3,2);
 
YMatrix_Port5 = ones(GroupsPerPort,1)*5;
 
 Vinf5 = Vinf(3,2);


%% port 6

GroupsPerPort = 4;
 Port6 = 0;
 ErrorPort6 = 0;
 Vinf6 = 0;
 
 Port6(1) = BL_thickness(3,3);
 ErrorPort6(1) = BL_error(3,3);
 
 Port6(2) = BL_thickness(4,3);
 ErrorPort6(2) = BL_error(4,3);
 
 Port6(3) = BL_thickness(1,4);
 ErrorPort6(3) = BL_error(1,4);
 
 Port6(4) = BL_thickness(2,4);
 ErrorPort6(4) = BL_error(2,4);
 
 
 Vinf6(1) = Vinf(3,3);
 Vinf6(2) = Vinf(4,3);
 Vinf6(3) = Vinf(1,3);
 Vinf6(4) = Vinf(2,3);
 
 
YMatrix_Port6 = ones(GroupsPerPort,1)*6;



%% port 7

GroupsPerPort = 4;
 Port7 = 0;
 ErrorPort7 = 0;
 
 Port7(1) = BL_thickness(5,1);
 ErrorPort7(1) = BL_error(5,1);
 
 Port7(2) = BL_thickness(6,1);
 ErrorPort7(2) = BL_error(6,1);
 
 Port7(3) = BL_thickness(7,3);
 ErrorPort7(3) = BL_error(7,3);
 
 Port7(4) = BL_thickness(8,3);
 ErrorPort7(4) = BL_error(8,3);
 
YMatrix_Port7 = ones(GroupsPerPort,1)*7;



 Vinf7 = Vinf(5,1);
 Vinf7 = Vinf(6,1);
 Vinf7 = Vinf(7,3);
 Vinf7 = Vinf(8,3);


%% port 8

GroupsPerPort = 2;
 Port8 = 0;
 ErrorPort8 = 0;
 
 Port8(1) = BL_thickness(5,2);
 ErrorPort8(1) = BL_error(5,2);
 
 Port8(2) = BL_thickness(6,2);
 ErrorPort8(2) = BL_error(6,2);

 Vinf8 = Vinf(5,2);
 Vinf8 = Vinf(6,2);

 
YMatrix_Port8 = ones(GroupsPerPort,1)*8;

%% port 9:

GroupsPerPort = 1;
Port9 = 0;
ErrorPort9 = 0;
 
Port9(1) = BL_thickness(5,3);
ErrorPort9(1) = BL_error(5,3);
 
 Vinf9 = Vinf(5,3);
 
YMatrix_Port9 = ones(GroupsPerPort,1)*9;


%% port 10

GroupsPerPort = 4;
 
Port10(1) = BL_thickness(7,1);
ErrorPort10(1) = BL_error(7,1);
 
Port10(2) = BL_thickness(8,1);
ErrorPort10(2) = BL_error(8,1);

Port10(3) = BL_thickness(5,4);
ErrorPort10(3) = BL_error(5,4);

Port10(4) = BL_thickness(6,4);
ErrorPort10(4) = BL_error(6,4);

 Vinf10 = Vinf(6,4);
 Vinf10 = Vinf(7,1);
 Vinf10 = Vinf(5,4);
 Vinf10 = Vinf(8,1);

YMatrix_Port10 = ones(GroupsPerPort,1)*10;


%% port 11


GroupsPerPort = 1;
 
Port11(1) = BL_thickness(8,2);
ErrorPort11(1) = BL_error(8,2);

Vinf11 = Vinf(8,2);

YMatrix_Port11 = ones(GroupsPerPort,1)*11;

%% get port location
PortLocation = [ 9.05 10.03 11.01 11.99 12.97 13.95 14.93 15.91 16.89 17.87 18.85 ] * 25.4 * 10^-3 ; %in m


%% Prediction 

% get free-stream conditions for our group

[ ignroe ignore2 ignore3 Vinfinity ] = Lab3Modified(11,1);

RynoldsFunc = @(x) ( (958.770603453 * Vinfinity ) / ( 1.8809518e-5 ) ) * x ;
TurblantFunc = @(x) ( 0.37*x ) / (RynoldsFunc(x))^(0.2); 
LaminarFunc = @(x) ( 5.2*x ) / sqrt((RynoldsFunc(x)));
%% all the ports: plots

%{ 

 turn on if you're using mean values

PortLocation = [ 9.05 10.03 11.01 11.99 12.97 13.95 14.93 15.91 16.89 17.87 18.85 ] %in inches
BLValues = [ Port1 Port2 Port3 Port4 WeightedAvg_p5 WeightedAvg_p6 WeightedAvg_p7 WeightedAvg_p8 WeightedAvg_p9 WeightedAvg_p10 WeightedAvg_p11 ]; %in mm
BLUncertinity = [ ErrorPort1 ErrorPort2 ErrorPort3 ErrorPort4 errorWeighted_p5 errorWeighted_p6 errorWeighted_p7 errorWeighted_p8 errorWeighted_p9 errorWeighted_p10 errorWeighted_p11 ];
%}


YMatrix_Port1 = YMatrix_Port1 .* PortLocation(1) ;
YMatrix_Port2 = YMatrix_Port2 .* PortLocation(2) ;
YMatrix_Port3 = YMatrix_Port3 .* PortLocation(3) ;
YMatrix_Port4 = YMatrix_Port4 .* PortLocation(4) ;
YMatrix_Port5 = YMatrix_Port5 .* PortLocation(5) ;
YMatrix_Port6 = YMatrix_Port6 .* PortLocation(6) ;
YMatrix_Port7 = YMatrix_Port7 .* PortLocation(7) ;
YMatrix_Port8 = YMatrix_Port8 .* PortLocation(8) ;
YMatrix_Port9 = YMatrix_Port9 .* PortLocation(9) ;
YMatrix_Port10 = YMatrix_Port10 .* PortLocation(10) ;
YMatrix_Port11 = YMatrix_Port11 .* PortLocation(11) ;


errorbar(YMatrix_Port1,Port1.*10^-3,ErrorPort1.*10^-3,'o')
hold on
errorbar(YMatrix_Port2,Port2.*10^-3,ErrorPort2.*10^-3,'o')
hold on
errorbar(YMatrix_Port3,Port3.*10^-3,ErrorPort3.*10^-3,'o')
hold on
errorbar(YMatrix_Port4,Port4.*10^-3,ErrorPort4.*10^-3,'o')
hold on
errorbar(YMatrix_Port5,Port5.*10^-3,ErrorPort5.*10^-3,'o')
hold on
errorbar(YMatrix_Port6,Port6.*10^-3,ErrorPort6.*10^-3,'o')
hold on
errorbar(YMatrix_Port7,Port7.*10^-3,ErrorPort7.*10^-3,'o')
hold on
errorbar(YMatrix_Port8,Port8.*10^-3,ErrorPort8.*10^-3,'o')
hold on
errorbar(YMatrix_Port9,Port9.*10^-3,ErrorPort9.*10^-3,'o')
hold on
errorbar(YMatrix_Port10,Port10.*10^-3,ErrorPort10.*10^-3,'o')
hold on
errorbar(YMatrix_Port11,Port11.*10^-3,ErrorPort11.*10^-3,'o')
hold on
fplot(TurblantFunc)
hold on
fplot(LaminarFunc)
hold off
xlim([-0.001 5.5])
ylim([-0.001 12e-3])
xlabel('X Location in m')
ylabel('Y location (Boundary layer thickness) in m')
title('Theoretical prediction Vs collected data for Boundary Layer')
legend('Port 1','Port 2','Port 3','Port 4','Port 5',...
'Port 6','Port 7','Port 8','Port 9','Port 10','Port 11',...
'Turbulent flow prediction','Laminar flow prediction');

grid minor


%% comparing the free-stream velocity as location of the probe changes
% to see how velocity changes


figure(2)
scatter( PortLocation , [ mean(Vinf1),mean(Vinf2),mean(Vinf3),mean(Vinf4),mean(Vinf5) , mean(Vinf6) ,mean(Vinf7),mean(Vinf8),mean(Vinf9),mean(Vinf10), mean(Vinf11) ] )
hold on
plot( PortLocation , [ mean(Vinf1),mean(Vinf2),mean(Vinf3),mean(Vinf4),mean(Vinf5) , mean(Vinf6) ,mean(Vinf7),mean(Vinf8),mean(Vinf9),mean(Vinf10), mean(Vinf11) ] )
title('Free-Stream Velocities with respect to port lcoatons')
xlabel('Port location (m)');
ylabel('Free-Stream Velocity (m/s)');
grid minor