%% Aero Lab 2 Main Function
%{

This script will compare groups and sections data for this lab, for more
info read Compare.m.


Authors:

    -Abdulla Alameri
    -Sam D'Souza
    -Foster Greer
    -Nathan Portman
    -Trevor Slack

%}
%% Housekeeping
clear;
clc;
close all;

%% Plot AOA versus Coefficient of lift (9m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,1)
        figure(1)
        hold on;
        leg = [ 'Section' num2str(j) ', Group ' num2str(p)]; %make a variable legend
        plot(alpha(1:3:end),Cl(1:3:end),'-*','DisplayName',leg)
        grid minor
        title(['\alpha vs C_L (Constant Airspeed 9 m/s)']);
        legend show
    end
end
xlabel('\alpha (Degrees)')
ylabel('C_L')
grid minor
hold off;
%% Plot AOA versus Coefficient of Drag (9m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,2)
        hold on;
        plot(alpha(1:3:end),Cd(1:3:end),'-*','DisplayName',leg) %plot every 3rd point for constant airspeed
        grid minor
        title(['\alpha vs C_D (Constant Airspeed 9 m/s)']);
    end
end
xlabel('\alpha (Degrees)')
ylabel('C_D')
grid minor
hold off;
%% Plot Coefficient of lift versus Coefficient of Drag (9m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,3)
        hold on;
        plot(Cl(1:3:end),Cd(1:3:end),'-*','DisplayName',leg)
        grid minor
        title(['C_L vs C_D (Constant Airspeed 9 m/s)']);
    end
end
xlabel('C_L')
ylabel('C_D')
grid minor
hold off;
%% Plot Cl/Cd verus AOA (9m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,4)
        hold on;
        plot(alpha(1:3:end),Cl(1:3:end)./Cd(1:3:end),'-*','DisplayName',leg)
        grid minor
        title(['\alpha vs C_L/C_D (Constant Airspeed 9 m/s)']);
    end
end
xlabel('\alpha (Degrees)')
ylabel('(C_L)/(C_D)')
grid minor
hold off;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot AOA versus Coefficient of lift (17m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        figure(2)
        subplot(2,2,1)
        hold on;
        leg = [ 'Section' num2str(j) ', Group ' num2str(p)]; %make a variable legend
        plot(alpha(2:3:end),Cl(2:3:end),'-*','DisplayName',leg)
        grid minor
        title(['\alpha vs C_L (Constant Airspeed 17 m/s)']);
        legend show
    end
end
xlabel('\alpha (Degrees)')
ylabel('C_L')
grid minor
hold off;
%% Plot AOA versus Coefficient of Drag (17m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,2)
        hold on;
        plot(alpha(2:3:end),Cd(2:3:end),'-*','DisplayName',leg) %plot every 3rd point for constant airspeed
        grid minor
        title(['\alpha vs C_D (Constant Airspeed 17 m/s)']);
    end
end
xlabel('\alpha (Degrees)')
ylabel('C_D')
grid minor
hold off;
%% Plot Coefficient of lift versus Coefficient of Drag (17m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,3)
        hold on;
        plot(Cl(2:3:end),Cd(2:3:end),'-*','DisplayName',leg)
        grid minor
        title(['C_L vs C_D (Constant Airspeed 17 m/s)']);
    end
end
xlabel('C_L')
ylabel('C_D')
grid minor
hold off;
%% Plot Cl/Cd verus AOA (17m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,4)
        hold on;
        plot(alpha(2:3:end),Cl(2:3:end)./Cd(2:3:end),'-*','DisplayName',leg)
        grid minor
        title(['\alpha vs C_L/C_D (Constant Airspeed 17 m/s)']);
    end
end
xlabel('\alpha (Degrees)')
ylabel('(C_L)/(C_D)')
grid minor
hold off;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot AOA versus Coefficient of lift (34m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        figure(3)
        subplot(2,2,1)
        hold on;
        leg = [ 'Section' num2str(j) ', Group ' num2str(p)]; %make a variable legend
        plot(alpha(3:3:end),Cl(3:3:end),'-*','DisplayName',leg)
        grid minor
        title(['\alpha vs C_L (Constant Airspeed 34 m/s)']);
        legend show
    end
end
xlabel('\alpha (Degrees)')
ylabel('C_L')
grid minor
hold off;
%% Plot AOA versus Coefficient of Drag (34m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,2)
        hold on;
        plot(alpha(3:3:end),Cd(3:3:end),'-*','DisplayName',leg) %plot every 3rd point for constant airspeed
        grid minor
        title(['\alpha vs C_D (Constant Airspeed 34 m/s)']);
    end
end
xlabel('\alpha (Degrees)')
ylabel('C_D')
grid minor
hold off;
%% Plot Coefficient of lift versus Coefficient of Drag (34m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,3)
        hold on;
        plot(Cl(3:3:end),Cd(3:3:end),'-*','DisplayName',leg)
        grid minor
        title(['C_L vs C_D (Constant Airspeed 34 m/s)']);
    end
end
xlabel('C_L')
ylabel('C_D')
grid minor
hold off;
%% Plot Cl/Cd verus AOA (34m/s)
for j=11:14
    for p=1:2:15
        [Cd,Cl,alpha]=Compare(j,p); %run windtunnel function
        subplot(2,2,4)
        hold on;
        plot(alpha(3:3:end),Cl(3:3:end)./Cd(3:3:end),'-*','DisplayName',leg)
        grid minor
        title(['\alpha vs C_L/C_D (Constant Airspeed 34 m/s)']);
    end
end
ylim([ -50 100])
xlabel('\alpha (Degrees)')
ylabel('(C_L)/(C_D)')
grid minor
hold off;