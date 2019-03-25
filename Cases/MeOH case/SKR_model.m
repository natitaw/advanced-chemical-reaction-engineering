function[fc_CH3OH, fc_CO, fc_CO2, fc_H2, fc_H2O]=SKR_model(T, P)
%==============================SRK model===================================
R= 8.314e-5; %Gas constant in m3*bar/K mol 
%---------------------------------------------Methanol Data---------------------------------------------------- 
%Data methanol 
T_ci_CH3OH = 512.6; %Critical temperature in Kelvin 
P_ci_CH3OH = 81; %Critical pressure in bar 
Acentric_CH3OH= 0.572; %Pitzer's acentric factor 
Tr_CH3OH= T./T_ci_CH3OH; %Reduced temperature 

%Calculations methanol constants  
m_methanol=0.48+1.574.*Acentric_CH3OH-0.176.*Acentric_CH3OH.^2;
alpha_methanol= (1+m_methanol.*(1-sqrt(Tr_CH3OH))).^2;
a_ci_methanol=0.42748.*alpha_methanol.*(R.*T_ci_CH3OH).^2./P_ci_CH3OH;
b_ci_methanol=0.08664.*R.*T_ci_CH3OH./P_ci_CH3OH;
A_methanol= (a_ci_methanol.*P)./(R.*T).^2;
B_methanol= (b_ci_methanol.*P)./(R.*T);

%---------------------------------------------Carbon monoxide calculations---------------------------------------------------- 
%Data Carbon monoxide
T_ci_CO = 132.9; %Critical temperature in Kelvin 
P_ci_CO = 35; %Critical pressure in bar 
Acentric_CO= 0.049; %Pitzer's acentric factor 
Tr_CO= T./T_ci_CO; %Reduced temperature 

%Calculations Carbon monoxide constants
m_CO=0.48+1.574.*Acentric_CO-0.176.*Acentric_CO.^2;
alpha_CO= (1+m_CO.*(1-sqrt(Tr_CO))).^2;
a_ci_CO=0.42748.*alpha_CO.*(R.*T_ci_CO).^2./P_ci_CO;
b_ci_CO=0.08664*R*T_ci_CO/P_ci_CO;
A_CO= (a_ci_CO*P)./(R.*T).^2;
B_CO= (b_ci_CO*P)./(R.*T);
%---------------------------------------------CO2 calculations---------------------------------------------------- 
%Data CO2 
T_ci_CO2 = 304.2; %Critical temperature in Kelvin 
P_ci_CO2 = 73.8; %Critical pressure in bar 
Acentric_CO2= 0.255; %Pitzer's acentric factor 
Tr_CO2= T./T_ci_CO2; %Reduced temperature 

%Calculations CO2 constants
m_CO2=0.48+1.574.*Acentric_CO2-0.176.*Acentric_CO2.^2;
alpha_CO2= (1+m_CO2.*(1-sqrt(Tr_CO2))).^2;
a_ci_CO2=0.42748.*alpha_CO2.*(R*T_ci_CO2).^2./P_ci_CO2;
b_ci_CO2=0.08664.*R.*T_ci_CO2./P_ci_CO2;
A_CO2= (a_ci_CO2.*P)./(R.*T)^2;
B_CO2= (b_ci_CO2.*P)./(R.*T);
%---------------------------------------------H2 calculations---------------------------------------------------- 
%Data H2 
T_ci_H2 = 43.6; %Critical temperature in Kelvin 
P_ci_H2 = 20.5; %Critical pressure in bar 
Acentric_H2= 0; %Pitzer's acentric factor 
Tr_H2= T./T_ci_H2; %Reduced temperature 

%Calculations H2 constants
m_H2=0.48+1.574.*Acentric_H2-0.176.*Acentric_H2.^2;
alpha_H2= (1+m_H2.*(1-sqrt(Tr_H2))).^2;
a_ci_H2=0.42748.*alpha_H2.*(R.*T_ci_H2).^2./P_ci_H2;
b_ci_H2=0.08664.*R.*T_ci_H2./P_ci_H2;
A_H2= (a_ci_H2.*P)./(R.*T).^2;
B_H2= (b_ci_H2.*P)./(R.*T);
%---------------------------------------------H2O calculations---------------------------------------------------- 
%Data H2O
T_ci_H2O = 647.3; %Critical temperature in Kelvin 
P_ci_H2O = 220.5; %Critical pressure in bar 
Acentric_H2O= 0.344; %Pitzer's acentric factor 
Tr_H2O= T./T_ci_H2O; %Reduced temperature 

%Calculations H2O constants
m_H2O=0.48+1.574.*Acentric_H2O-0.176.*Acentric_H2O.^2;
alpha_H2O= (1+m_H2O.*(1-sqrt(Tr_H2O))).^2;
a_ci_H2O=0.42748.*alpha_H2O.*(R.*T_ci_H2O).^2./P_ci_H2O;
b_ci_H2O=0.08664.*R.*T_ci_H2O./P_ci_H2O;
A_H2O= (a_ci_H2O.*P)./(R.*T).^2;
B_H2O= (b_ci_H2O.*P)./(R.*T);  

%-----------Z-Factor calculation---------------------------------------------------------- 
%Methanol
Zv_methanol_all=(roots([1 -1 A_methanol-B_methanol-B_methanol.^2 -A_methanol.*B_methanol]));
if isreal(Zv_methanol_all(1))~=0
    b=1;
elseif isreal(Zv_methanol_all(2))~=0
    b=2;
elseif isreal(Zv_methanol_all(2))==0
    b=3;
end
Zv_methanol_V=Zv_methanol_all(b);
log_f_p_CH3OH=Zv_methanol_V-1-log(Zv_methanol_V-B_methanol)-(A_methanol/B_methanol).*log((Zv_methanol_V+B_methanol)./Zv_methanol_V);
%CO
Zv_CO_all=(roots([1 -1 A_CO-B_CO-B_CO.^2 -A_CO.*B_CO]));
if isreal(Zv_CO_all(1))~=0
    b2=1;
elseif isreal(Zv_CO_all(2))~=0
    b2=2;
elseif isreal(Zv_CO_all(2))==0
    b2=3;
end
Zv_CO_V=Zv_CO_all(b2);
log_f_p_CO=Zv_CO_V-1-log(Zv_CO_V-B_CO)-(A_CO./B_CO).*log((Zv_CO_V+B_CO)./Zv_CO_V);
%CO2
Zv_CO2_all=(roots([1 -1 A_CO2-B_CO2-B_CO2.^2 -A_CO2.*B_CO2]));
if isreal(Zv_CO2_all(1))~=0
    b3=1;
elseif isreal(Zv_CO2_all(2))~=0
    b3=2;
elseif isreal(Zv_CO2_all(2))==0
    b3=3;
end
Zv_CO2_V=Zv_CO2_all(b3);
log_f_p_CO2=Zv_CO2_V-1-log(Zv_CO2_V-B_CO2)-(A_CO2./B_CO2).*log((Zv_CO2_V+B_CO2)./Zv_CO2_V);
%H2
Zv_H2_all=(roots([1 -1 A_H2-B_H2-B_H2.^2 -A_H2.*B_H2]));
if isreal(Zv_H2_all(1))~=0
    b4=1;
elseif isreal(Zv_H2_all(2))~=0
    b4=2;    
elseif isreal(Zv_H2_all(2))==0
    b4=3;
end
Zv_H2_V=Zv_H2_all(b4);
log_f_p_H2=Zv_H2_V-1-log(Zv_H2_V-B_H2)-(A_H2/B_H2).*log((Zv_H2_V+B_H2)./Zv_H2_V);
%H2O
Zv_H2O_all=(roots([1 -1 A_H2O-B_H2O-B_H2O.^2 -A_H2O.*B_H2O]));
if isreal(Zv_H2O_all(1))~=0
    b1=1;
elseif isreal(Zv_H2O_all(1))~=0
    b1=2;    
elseif isreal(Zv_H2O_all(2))==0
    b1=3;
end
Zv_H2O_V=Zv_H2O_all(b1);
log_f_p_H2O=Zv_H2O_V-1-log(Zv_H2O_V-B_H2O)-(A_H2O/B_H2O).*log((Zv_H2O_V+B_H2O)./Zv_H2O_V);
%===================Fugacity coefficients==============================
fc_CH3OH=exp(log_f_p_CH3OH);
fc_CO=exp(log_f_p_CO);
fc_CO2=exp(log_f_p_CO2);
fc_H2=exp(log_f_p_H2);
fc_H2O=exp(log_f_p_H2O);
end