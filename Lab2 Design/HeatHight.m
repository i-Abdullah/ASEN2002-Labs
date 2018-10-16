function [h] = HeatHight(heat,l,u)
% ---------------------------------------------------------------
% October 11th, 2018
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

Interval = l:0.01:u;

[ Temp Sound Press DensityValues ] = atmoscoesa(Interval) ;

low = 1;
high = length(DensityValues);
mid = floor((low+high)/2);
tol = 1e-5; %tolerance kg/m^3

while abs(DensityValues(high) - DensityValues(low)) > tol
    
    if DensityValues(mid) > heat
        
        low = mid;
        high = high;
        mid = floor((low+high)/2);
        
        
        
    else
        
        high = mid;
        low = low;
        mid = floor((low+high)/2);

    end
    
    
end

h = floor((low+high)/2)/100 ;

end