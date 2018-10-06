%% Ideal Gas Law

%{
P = ;
V = ;
R = 8.314; % J/ (mol * k)
T = ;

n = (P*V)/(R*T);

%}
%% Uncertinty:

m_maylar = 9.80
uncert_maylar = 0.005;
m_He = 1.9304622;
uncert_He = 0.004;
m_load = 0.525;
uncert_load = 0.006;

Total_mass = m_maylar + m_He + m_load;
Total_mass_uncert = uncert_He + uncert_load + uncert_maylar;

% Fractional Mass:
massFract_He = m_He / Total_mass ;
massFract_He_unct = sqrt((uncert_He/m_He)^2+ (Total_mass_uncert/Total_mass)^2) * massFract_He;

massFract_load = m_load/Total_mass;
massFract_load_unct = sqrt((uncert_load/m_load)^2+ (Total_mass_uncert/Total_mass)^2) * massFract_load;

massFract_maylar = m_maylar/Total_mass;
massFract_maylar_unct = sqrt((uncert_maylar/m_maylar)^2+ (Total_mass_uncert/Total_mass)^2) * massFract_maylar;