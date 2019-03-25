function[rCH3OH_A3, rH2O_B2, rCH3OH_C3]=methanol_kinetics_gas(yCH3OH, yCO, yCO2, yH2,yH2O, T, P, dp)
%Kinetics of the methanol synthesis reactions
%Thiele_model, SKR_model and Equilibrium_model required in path to run
%Required inputs are the gas fractions, temperature (K), pressure (bar) and
%particle diameter in (m) 
%This file calculates the rates of the methanol synthesis reactions in mol/kg_cat*s  
%Rates include the intra-particle diffusion limitations 

R=8.314; %Gas constant in J K-1 mol
%===============================-kinetics Gas=================================
K00p1= 2.391e-13; 
K00p2= 1.068e2; 
K00p3= 2.554e-11;
k0_psA3= 4.89e7;
k0_psB2= 9.64e11;
k0_psC3= 1.09e5;
K0CO= 2.16e-5;
k0CO2= 7.05e-7;
K0H2O_KH2= 6.37e-9; 
%===========================Activation energies============================ 
E_K0p1= 98388; 
E_K0p2= -39683;  
E_K0p3= 58705; 
E_k0_psA3= -113000;
E_k0_psB2= -152900;
E_k0_psC3= -87500; 
E_K0CO= 46800; 
E_k0CO2= 61700; 
E_K0H2O_KH2= 84000; 

k0p1= K00p1*exp(E_K0p1/(R*T));
k0p2= K00p2*exp(E_K0p2/(R*T));
k0p3= K00p3*exp(E_K0p3/(R*T));
k_psA3= k0_psA3*exp(E_k0_psA3/(R*T));
k_psB2= k0_psB2*exp(E_k0_psB2/(R*T));
k_psC3= k0_psC3*exp(E_k0_psC3/(R*T));
KCO= K0CO*exp(E_K0CO/(R*T));
KCO2= k0CO2*exp(E_k0CO2/(R*T));
KH2O_KH2= K0H2O_KH2*exp(E_K0H2O_KH2/(R*T));

%===============================Fugacities (bar) ===============================
[fc_CH3OH, fc_CO, fc_CO2, fc_H2, fc_H2O]=SKR_model(T, P); 
f_CH3OH=fc_CH3OH.*yCH3OH.*P;
f_CO=fc_CO.*yCO.*P;
f_CO2=fc_CO2.*yCO2.*P;
f_H2=fc_H2.*yH2.*P;
f_H2O=fc_H2O.*yH2O.*P;

%=========================Rates (%mol/m3*s) two phase without intra particle difffusion limitations================================
Demon=((1+KCO*f_CO+KCO2*f_CO2)*(sqrt(f_H2)+(KH2O_KH2)*f_H2O));
rCH3OH_A3_bulk= ((k_psA3*KCO*(f_CO*f_H2^1.5-f_CH3OH/(sqrt(f_H2)*k0p1)))/Demon); 
rH2O_B2_bulk=  ((k_psB2*KCO2*(f_CO2*f_H2-f_H2O*f_CO/k0p2))/Demon);
rCH3OH_C3_bulk= ((k_psC3*KCO2*(f_CO2*f_H2^(1.5)-f_CH3OH*f_H2O/((f_H2^(3/2))*k0p3)))/Demon);

%=======================Effectiveness factor===============================
[yCOeq, yH2eq, yCO2eq, yH2Oeq, yCH3OHeq]=Equilibrium_model(yCH3OH, yCO, yCO2, yH2,yH2O, T, P);
Keq_pseudo=(yCH3OHeq/yH2Oeq);

[Eff_factor]=Thiele_model(Keq_pseudo, rCH3OH_A3_bulk, rH2O_B2_bulk, rCH3OH_C3_bulk,yCH3OH, yCO, yCO2, yH2,yH2O, T, P, dp);

%======================Rates two phass including diffusion limitations=====
rCH3OH_A3=rCH3OH_A3_bulk.*Eff_factor; %mol/kg_cat*s
rH2O_B2=rH2O_B2_bulk.*Eff_factor; %mol/kg_cat*s
rCH3OH_C3=rCH3OH_C3_bulk.*Eff_factor; %mol/kg_cat*s
end
