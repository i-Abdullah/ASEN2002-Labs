function [] = modelAir(hight,err)
%% Intro
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
% The Following function is going to model and plot the Atmospheric
% conditions at a certain hight and around a bound of error. 
%
%
% ---------------------------------------------------------------

%Esablish a range of values.

h = ((hight-(err+100)):(hight+err+100));

[T,a,P,rho] = atmoscoesa(h,'Error');
%error will through error if hight exceeds the hight the atm model can predict.

%% Plot the Results in subplot
subplot(2,2,[1 2])

ylim ([min(h) max(h)])
plot(T,h,'-.b', 'LineWidth',1.5)
xlabel('Temperature (K)');
ylabel('Altitude (m)');
title('Temperature Changes with Altitude');
grid minor
hold on
ref = refline(0,hight)
ref.Color = 'r';
hold on
ref = refline(0,min(h)+100)
ref.Color = 'r';
hold on
ref = refline(0,max(h)-100)
ref.Color = 'r';
ylim ([min(h) max(h)])
yticks ([min(h):200:max(h)])
xticks ([min(T):2:max(T)])
legend('h VS T','Upper error','Hight Goal','Lower error')
hold off


subplot(2,2,3)
plot(P./1000,h,'-.b', 'LineWidth',1.5)
xlabel('Pressure (kPa)');
ylabel('Altitude (m)');
title('Pressure Changes with Altitude');
grid minor
ylim ([min(h) max(h)])
hold on
ref = refline(0,hight)
ref.Color = 'r';
hold on
ref = refline(0,min(h)+100)
ref.Color = 'r';
hold on
ref = refline(0,max(h)-100)
ref.Color = 'r';
ylim ([min(h) max(h)])
yticks ([min(h):200:max(h)])
legend('h VS P','Upper error','Hight Goal','Lower error')
hold off

subplot(2,2,4)
plot(rho,h,'-.b', 'LineWidth',1.5)
xlabel('Density (Kg/m^3)');
ylabel('Altitude (m)');
title('Density Changes with Altitude');

grid minor
ylim ([min(h) max(h)])
hold on
ref = refline(0,hight)
ref.Color = 'r';
hold on
ref = refline(0,min(h)+100)
ref.Color = 'r';
hold on
ref = refline(0,max(h)-100)
ref.Color = 'r';
ylim ([min(h) max(h)])
yticks ([min(h):200:max(h)])
legend('h VS rho','Upper error','Hight Goal','Lower error')
hold off

end