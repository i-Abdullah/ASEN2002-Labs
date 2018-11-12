function [ sigmaVPitot] = SigmaV(R,Airdiff,Tatm,Patm,x,sigmaD_P,sigma_P,sigma_T)
%% Calculate the uncertainty in the velocity in the Pitot tube 

% Create a symbolic function
syms Air_diff T_atm P_atm

AirstreamPitot(Air_diff,T_atm,P_atm) = sqrt((2*Air_diff*R*T_atm)/(P_atm)); 


AirstreamPitotx = zeros(1,x(1));
VuncertaintyPitot = zeros(1,x(1));

for i = 1:x(1)
   AirstreamPitotx(i) =  double(AirstreamPitot(Airdiff(i),Tatm(i),Patm(i))); 
end

DDPPitot = diff(AirstreamPitot,Air_diff);
DTPitot= diff(AirstreamPitot,T_atm);
DPPitot= diff(AirstreamPitot,P_atm);

for i = 1:x(1)
VuncertaintyPitot(i) = sqrt((double(DDPPitot(Airdiff(i),Tatm(i),Patm(i)))*sigmaD_P(i))^2+...
    (double(DTPitot(Airdiff(i),Tatm(i),Patm(i)))*sigma_P(i))^2+...
    (double(DPPitot(Airdiff(i),Tatm(i),Patm(i)))*sigma_T(i))^2);
end

sigmaVPitot = VuncertaintyPitot;

%{
% Uncertainty in the atmospheric pressure for each voltage in Pa 
sigmaP1 = 0.01*patm1;
sigmaP2 = 0.01*patm2;
sigmaP3 = 0.01*patm3;
sigmaP4 = 0.01*patm4;
sigmaP5 = 0.01*patm5;

% Uncertainty in the differential pressure in Pa
sigmaDP = 0.01*.24884; 

% Uncertainty in the atmospheric temperature in K
sigmaT = 1/4;

% Partial derivative of dvddP (Manometer)
dvddP1Mano = (R.*temp1)./(sqrt(2).*(1-A1A2ratio^2).*patm1.*sqrt(R.*temp1.*ManoDiffPressue1./((1-A1A2ratio^2).*patm1)));
dvddP2Mano = (R.*temp2)./(sqrt(2).*(1-A1A2ratio^2).*patm2.*sqrt(R.*temp2.*ManoDiffPressue2./((1-A1A2ratio^2).*patm2)));
dvddP3Mano = (R.*temp3)./(sqrt(2).*(1-A1A2ratio^2).*patm3.*sqrt(R.*temp3.*ManoDiffPressue3./((1-A1A2ratio^2).*patm3)));
dvddP4Mano = (R.*temp4)./(sqrt(2).*(1-A1A2ratio^2).*patm4.*sqrt(R.*temp4.*ManoDiffPressue4./((1-A1A2ratio^2).*patm4)));
dvddP5Mano = (R.*temp5)./(sqrt(2).*(1-A1A2ratio^2).*patm5.*sqrt(R.*temp5.*ManoDiffPressue5./((1-A1A2ratio^2).*patm5)));

% Partial derivative of dvddP (ASD)
dvddP1ASD = (R.*temp1)./(sqrt(2).*(1-A1A2ratio^2).*patm1.*sqrt(R.*temp1.*diffPressure1./((1-A1A2ratio^2).*patm1)));
dvddP2ASD = (R.*temp2)./(sqrt(2).*(1-A1A2ratio^2).*patm2.*sqrt(R.*temp2.*diffPressure2./((1-A1A2ratio^2).*patm2)));
dvddP3ASD = (R.*temp3)./(sqrt(2).*(1-A1A2ratio^2).*patm3.*sqrt(R.*temp3.*diffPressure3./((1-A1A2ratio^2).*patm3)));
dvddP4ASD = (R.*temp4)./(sqrt(2).*(1-A1A2ratio^2).*patm4.*sqrt(R.*temp4.*diffPressure4./((1-A1A2ratio^2).*patm4)));
dvddP5ASD = (R.*temp5)./(sqrt(2).*(1-A1A2ratio^2).*patm5.*sqrt(R.*temp5.*diffPressure5./((1-A1A2ratio^2).*patm5)));


% Partial derivative of dvdT (Manometer)
dvdT1Mano = (R.*ManoDiffPressue1)./(sqrt(2).*(1-A1A2ratio^2).*patm1.*sqrt(R.*temp1.*ManoDiffPressue1./((1-A1A2ratio^2).*patm1)));
dvdT2Mano = (R.*ManoDiffPressue2)./(sqrt(2).*(1-A1A2ratio^2).*patm2.*sqrt(R.*temp2.*ManoDiffPressue2./((1-A1A2ratio^2).*patm2)));
dvdT3Mano = (R.*ManoDiffPressue3)./(sqrt(2).*(1-A1A2ratio^2).*patm3.*sqrt(R.*temp3.*ManoDiffPressue3./((1-A1A2ratio^2).*patm3)));
dvdT4Mano = (R.*ManoDiffPressue4)./(sqrt(2).*(1-A1A2ratio^2).*patm4.*sqrt(R.*temp4.*ManoDiffPressue4./((1-A1A2ratio^2).*patm4)));
dvdT5Mano = (R.*ManoDiffPressue5)./(sqrt(2).*(1-A1A2ratio^2).*patm5.*sqrt(R.*temp5.*ManoDiffPressue5./((1-A1A2ratio^2).*patm5)));

% Partial derivative of dvdT (ASD)
dvdT1ASD = (R.*diffPressure1)./(sqrt(2).*(1-A1A2ratio^2).*patm1.*sqrt(R.*temp1.*diffPressure1./((1-A1A2ratio^2).*patm1)));
dvdT2ASD = (R.*diffPressure2)./(sqrt(2).*(1-A1A2ratio^2).*patm2.*sqrt(R.*temp2.*diffPressure2./((1-A1A2ratio^2).*patm2)));
dvdT3ASD = (R.*diffPressure3)./(sqrt(2).*(1-A1A2ratio^2).*patm3.*sqrt(R.*temp3.*diffPressure3./((1-A1A2ratio^2).*patm3)));
dvdT4ASD = (R.*diffPressure4)./(sqrt(2).*(1-A1A2ratio^2).*patm4.*sqrt(R.*temp4.*diffPressure4./((1-A1A2ratio^2).*patm4)));
dvdT5ASD = (R.*diffPressure5)./(sqrt(2).*(1-A1A2ratio^2).*patm5.*sqrt(R.*temp5.*diffPressure5./((1-A1A2ratio^2).*patm5)));

% Parial derivative of dvdP (Manometer)
dvdP1Mano = (-ManoDiffPressue1.*R.*temp1)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp1.*ManoDiffPressue1./((1-A1A2ratio^2).*patm1)).*patm1.^2);
dvdP2Mano = (-ManoDiffPressue2.*R.*temp2)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp2.*ManoDiffPressue2./((1-A1A2ratio^2).*patm2)).*patm2.^2);
dvdP3Mano = (-ManoDiffPressue3.*R.*temp3)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp3.*ManoDiffPressue3./((1-A1A2ratio^2).*patm3)).*patm3.^2);
dvdP4Mano = (-ManoDiffPressue4.*R.*temp4)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp4.*ManoDiffPressue4./((1-A1A2ratio^2).*patm4)).*patm4.^2);
dvdP5Mano = (-ManoDiffPressue5.*R.*temp5)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp5.*ManoDiffPressue5./((1-A1A2ratio^2).*patm5)).*patm5.^2);

% Parial derivative of dvdP (ASD)
dvdP1ASD = (-diffPressure1.*R.*temp1)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp1.*diffPressure1./((1-A1A2ratio^2).*patm1)).*patm1.^2);
dvdP2ASD = (-diffPressure2.*R.*temp2)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp2.*diffPressure2./((1-A1A2ratio^2).*patm2)).*patm2.^2);
dvdP3ASD = (-diffPressure3.*R.*temp3)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp3.*diffPressure3./((1-A1A2ratio^2).*patm3)).*patm3.^2);
dvdP4ASD = (-diffPressure4.*R.*temp4)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp4.*diffPressure4./((1-A1A2ratio^2).*patm4)).*patm4.^2);
dvdP5ASD = (-diffPressure5.*R.*temp5)./(sqrt(2).*(1-A1A2ratio^2).*sqrt(R.*temp5.*diffPressure5./((1-A1A2ratio^2).*patm5)).*patm5.^2);


%uncertainty in velocity (Manometer)
sigmav1Mano = sqrt((dvddP1Mano.*sigmaDP).^2+(dvdT1Mano.*sigmaT).^2+(dvdP1Mano.*sigmaP1).^2);
sigmav2Mano = sqrt((dvddP2Mano.*sigmaDP).^2+(dvdT2Mano.*sigmaT).^2+(dvdP2Mano.*sigmaP2).^2);
sigmav3Mano = sqrt((dvddP3Mano.*sigmaDP).^2+(dvdT3Mano.*sigmaT).^2+(dvdP3Mano.*sigmaP3).^2);
sigmav4Mano = sqrt((dvddP4Mano.*sigmaDP).^2+(dvdT4Mano.*sigmaT).^2+(dvdP4Mano.*sigmaP4).^2);
sigmav5Mano = sqrt((dvddP5Mano.*sigmaDP).^2+(dvdT5Mano.*sigmaT).^2+(dvdP5Mano.*sigmaP5).^2);

%uncertainty in velocity (ASD)
sigmav1ASD = sqrt((dvddP1ASD.*sigmaDP).^2+(dvdT1ASD.*sigmaT).^2+(dvdP1ASD.*sigmaP1).^2);
sigmav2ASD = sqrt((dvddP2ASD.*sigmaDP).^2+(dvdT2ASD.*sigmaT).^2+(dvdP2ASD.*sigmaP2).^2);
sigmav3ASD = sqrt((dvddP3ASD.*sigmaDP).^2+(dvdT3ASD.*sigmaT).^2+(dvdP3ASD.*sigmaP3).^2);
sigmav4ASD = sqrt((dvddP4ASD.*sigmaDP).^2+(dvdT4ASD.*sigmaT).^2+(dvdP4ASD.*sigmaP4).^2);
sigmav5ASD = sqrt((dvddP5ASD.*sigmaDP).^2+(dvdT5ASD.*sigmaT).^2+(dvdP5ASD.*sigmaP5).^2);

% Outputs 
sigmaVMano = [sigmav1Mano,sigmav2Mano,sigmav3Mano,sigmav4Mano,sigmav5Mano];
sigmaVASD = [sigmav1ASD,sigmav2ASD,sigmav3ASD,sigmav4ASD,sigmav5ASD];
    
    %}
end
