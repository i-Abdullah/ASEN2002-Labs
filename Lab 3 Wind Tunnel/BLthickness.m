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
            
       [ ignored, Thickness, Error ]  = Lab3Modified(i,groups(j));
        BL_thickness(j,i-10) = Thickness;
        BL_error(j,i-10) = Error;

        end            
  
    
    end
    
end

% there was total 11 ports, more than one group connected to each port.

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
 
for i=5:6
   Weight = ( 1/(BL_error(i,2)^2) );
   WeightedAvg_p9 = WeightedAvg_p9 + ( (Weight*BL_thickness(i,2))) ;
   errorWeighted_p9 = errorWeighted_p9 + Weight;
end



WeightedAvg_p9 = WeightedAvg_p9 /errorWeighted_p9;
errorWeighted_p9 = 2 / sqrt((errorWeighted_p9));
