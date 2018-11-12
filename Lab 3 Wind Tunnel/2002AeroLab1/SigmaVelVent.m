function [ sigmaVVent ] = SigmaVelVent(R,A1A2ratio,ManoDiffPressueTotal,Tatm,Patm,x,sigmaD_P,sigma_P,sigma_T)

syms Mano_DiffPressue Tatmosphere Patmosphere

AirstreamVent(Mano_DiffPressue,Tatmosphere,Patmosphere) = sqrt( (2 * Mano_DiffPressue * R * Tatmosphere) / ( (Patmosphere) * (1-(A1A2ratio)^2) ) );

VuncertaintyVent = zeros(1,x(1));

for i = 1:x(1)
   AirstreamVentx(i) =  double(AirstreamVent(ManoDiffPressueTotal(i),Tatm(i),Patm(i))) ; 
end

DDPVent = diff(AirstreamVent,Mano_DiffPressue);
DTVent = diff(AirstreamVent,Tatmosphere);
DPVent = diff(AirstreamVent,Patmosphere);

for i = 1:x(1)
VuncertaintyVent(i) = sqrt((double(DDPVent(ManoDiffPressueTotal(i),Tatm(i),Patm(i)))*sigmaD_P(i))^2+...
    (double(DTVent(ManoDiffPressueTotal(i),Tatm(i),Patm(i)))*sigma_P(i))^2+...
    (double(DPVent(ManoDiffPressueTotal(i),Tatm(i),Patm(i)))*sigma_T(i))^2);
end

sigmaVVent  = VuncertaintyVent;
end

