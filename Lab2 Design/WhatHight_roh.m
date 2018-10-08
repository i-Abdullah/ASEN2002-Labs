function h = WhatHight_roh(roh)
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
%           - Density of your fluid (kg/m^3)
%
% ---------------------------------------------------------------
% OUTPUTS:
%           - Hight in m.

roh_air = 10;
i = 0;

while abs(roh_air - roh) >= 1e-4
    %the more you want to be accurate with hight you'll change the
    %conditioon, however your input hight must be more accurate too!
    
    %Unless the density is inputted to a very high degree of accuracy
    %relative to how many sig figs, this function need some bigger
    %difference. 7 is arbitrary chosen but it is a good choice.
    
    i = i+1; %start by 1 meter, and keep adding a meter everytime just to get our measurments accurate enough.
   
    [ a b c d ] = atmosisa(i); %we get our data @ each hight.
    
    roh_air = d; %the output named c corresponds to the density of the air.
    
end

h = i; %h, our hight will be basically how many i's we have, each i represent 1 m.

end