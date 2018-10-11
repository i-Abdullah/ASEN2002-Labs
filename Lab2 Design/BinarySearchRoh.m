function [h] = BinarySearchRoh(densityGiven,l,u)
% ---------------------------------------------------------------
% October 5th, 2018
% Done by :
%              1- Mia Abouhamad
%              2- Abdulla Alameri
%              3- Daniel Barth
%              4- Zihao Ding
%              5- Chava Friedman
%              6- Eric Hunnel
%              7- Vinay Simlot
%              8- Ryan Smithers
%
% ---------------------------------------------------------------
%
% This function is designed to figure out a hight utilizing the atmospheric
% model. This function will do the math backward from a given density to
% find the hight you should be at if your density changed holding
% everything else to be constant. We know that the density inside our
% ballon will be equal to the density of the air outside the ballon when
% the ballon is in neutral buoyancy. Thus, the function will go and try to
% match the same density and tell us what hight is that.
%
% ---------------------------------------------------------------
% INPUTS:
%           - densityGiven : target density
%           - l : Lower bound of hight that you think the density lies within
%           - u : upper bound of hight you think the density lies within
%
% ---------------------------------------------------------------
% OUTPUTS:
%           - Hight in m.

Interval = l:0.1:u;

[ Temp Sound Press DensityValues ] = atmoscoesa(Interval) ;

HalfIndex = (length(DensityValues)/2);
halfWay = DensityValues(HalfIndex);


%Since density keeps increasing as we go up, our condition will be 


while abs (DensityValues(1) - DensityValues(end)) > 1e-3
    if(length(DensityValues)<3)
        break
    end
    if halfWay > densityGiven
              
        l = l;
        u = halfWay;
        DensityValues = DensityValues(1:halfWay);
        halfWay = DensityValues((length(DensityValues)/2));
        HalfIndex = (length(DensityValues)/2);
        halfWay = DensityValues(HalfIndex);


    else
        
        u = u;
        l = halfWay;
       DensityValues = DensityValues(index:end);
        Interval = 1:length(DensityValues);




        
    end
    %{
 Interval = l:0.1:u;
 
 
 if Interval(1) == 0
DensityValues = DensityValues(1:length(Interval));

 else
     
     DensityValues = DensityValues(Interval(1):length(Interval));

     
 end
    %}

halfWay = DensityValues(length(DensityValues)/2);

    
end

h = (l+u)/2

end