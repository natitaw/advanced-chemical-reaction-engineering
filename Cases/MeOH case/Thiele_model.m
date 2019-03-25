function[Eff_factor]=Thiele_model(Keq_pseudo, rCH3OH_A3_bulk, rH2O_B2_bulk, rCH3OH_C3_bulk, yCH3OH, yCO, yCO2, yH2,yH2O, T, P, dp)
R=8.314; %Gas constant in J K-1 mol
R_bar= 8.314e-5; %Gas constant in m3*bar/K mol 
Ctot=P/(R_bar*T);
rp=dp./2; %Radius particle in m
rho_cat= 1950; %kg/m3 

%----------Diffusion calculations--------------------------
a= 10e-9; % mean pore radius 
epsilon_tau= 0.123; %ratio porosity 
Patm=P*0.986923267;
M_methanol= 32.04; %Molar weightin g/mol
M_H2O= 18.02; %Molar weightin g/mol
M_H2= 2.016; %Molar weightin g/mol
M_CO= 28.01; %Molar weightin g/mol
M_CO2= 44.01; %Molar weightin g/mol

Deff_methanol=a*epsilon_tau*(2/3)*sqrt((8*R*T)/(pi*M_methanol/1000));

%---------------Diffusion volumes of simple molecules---------------
Dv_methanol=29.9; %cm3/mol
Dv_H2=7.07;%cm3/mol
Dv_H2O=12.7;%cm3/mol
Dv_CO=18.9;%cm3/mol
Dv_CO2=26.9;%cm3/mol

D_methanol_H2=epsilon_tau*((1e-3.*T.^(1.75).*(M_methanol.^-1+M_H2.^-1).^0.5)./(Patm.*(Dv_methanol.^(1./3)+Dv_H2.^(1./3)).^2))./100;
D_methanol_H2O=epsilon_tau*((1e-3.*T.^(1.75).*(M_methanol.^-1+M_H2O.^-1).^0.5)./(Patm.*(Dv_methanol.^(1./3)+Dv_H2O.^(1./3)).^2))./100;
D_methanol_CO=epsilon_tau*((1e-3.*T.^(1.75).*(M_methanol.^-1+M_CO.^-1).^0.5)./(Patm.*(Dv_methanol.^(1./3)+Dv_CO.^(1./3)).^2))./100;
D_methanol_CO2=epsilon_tau*((1e-3.*T.^(1.75).*(M_methanol.^-1+M_CO2.^-1).^0.5)./(Patm.*(Dv_methanol.^(1./3)+Dv_CO2.^(1./3)).^2))./100;

%----------------Realistic diffusion---------------------------------------
Dcal_methanol=1/Deff_methanol+yH2/D_methanol_H2+yH2O/D_methanol_H2O+yCO/D_methanol_CO+yCO2/D_methanol_CO2;

Drel_methanol=1/Dcal_methanol;

r_methanol=(rCH3OH_A3_bulk.*rho_cat+rCH3OH_C3_bulk.*rho_cat);

Driving_force=(Ctot*yH2-Ctot*yCH3OH/Keq_pseudo);

ki_pseudo=r_methanol/Driving_force;
       
phi_m=(rp./2.5).*sqrt((ki_pseudo.*(Keq_pseudo+1))./(Drel_methanol.*Keq_pseudo));
Eff_factor=(1./phi_m).*(((3.*phi_m).*coth(3.*phi_m)-1)./(3.*phi_m));
end