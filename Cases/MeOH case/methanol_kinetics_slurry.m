function[rCH3OH_A3, rH2O_B2, rCH3OH_C3]=methanol_kinetics_slurry(C_CH3OH, C_CO, C_CO2, C_H2,C_H2O, T)
%Kinetics of the methanol synthesis reactions
%Required inputs are the concentrations in the slurry phase in mol/m3, temperature (K)
%This file calculates the rates of the methanol synthesis reactions in mol/kg_cat*s

R=8.314; %Gas constant in J K-1 mol
%===============================-kinetics Gas=================================
K00p1= 1.72e-16; 
K00p2= 5.81e1; 
K00p3= 9.99e-15;
k0_psA3=  3.824e6;
k0_psB2= 2.696e11;
k0_psC3= 8.524e3;
K0CO= 3.78e-6;
k0CO2= 2.834e-7;
K0H2O_KH2= 7.52e-9; 
%===========================Activation energies============================ 
E_K0p1= 126011; 
E_K0p2= -33760;  
E_K0p3= 92251; 
E_k0_psA3= -108125;
E_k0_psB2= -150463;
E_k0_psC3= -82625; 
E_K0CO= 47438; 
E_k0CO2= 54753; 
E_K0H2O_KH2= 72930; 

k0p1= K00p1*exp(E_K0p1/(R*T));
k0p2= K00p2*exp(E_K0p2/(R*T));
k0p3= K00p3*exp(E_K0p3/(R*T));
k_psA3= k0_psA3*exp(E_k0_psA3/(R*T));
k_psB2= k0_psB2*exp(E_k0_psB2/(R*T));
k_psC3= k0_psC3*exp(E_k0_psC3/(R*T));
KCO= K0CO*exp(E_K0CO/(R*T));
KCO2= k0CO2*exp(E_k0CO2/(R*T));
KH2O_KH2= K0H2O_KH2*exp(E_K0H2O_KH2/(R*T));

%=========================Rates two phase================================
Demon=((1+KCO*C_CO+KCO2*C_CO2)*(sqrt(C_H2)+(KH2O_KH2)*C_H2O));
rCH3OH_A3= ((k_psA3*KCO*(C_CO*C_H2^1.5-C_CH3OH/(sqrt(C_H2)*k0p1)))/Demon); %mol/m3*s
rH2O_B2=  ((k_psB2*KCO2*(C_CO2*C_H2-C_H2O*C_CO/k0p2))/Demon);
rCH3OH_C3= ((k_psC3*KCO2*(C_CO2*C_H2^(1.5)-C_CH3OH*C_H2O/((C_H2^(3/2))*k0p3)))/Demon);
 
end
