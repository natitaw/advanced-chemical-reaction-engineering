 function []=eo_reactor1D
%==========================================================================
%2016 EO reactor, S.R.A. Kersten
%1D model
%assumptions: velocity is constant, no pressure drop
%==========================================================================
   



%the variables below can be changed in the simulations----
%--------------variables----------------------------------
u=1.3;      %m/s - superficial gas velocity == constant
L=30;       %m reactor length
Tin=225;    %C - inlet temperature
P=1e6;      %Pa - pressure == constant 
fo2=0.12;   %vol. fraction - vol. f. O2 in the feed
fe=1-fo2;   %vol. fraction - vol. f. E in the feed - the feed contains only O2 and E
dp=4e-3;    %m - diameter of the catalyst particles
dtube=2 ;   %inch - tube diameter, only 1, 2 and 30 inch can be selected
%----------------------------------------------------------

%--------------fixed paramters---------------------------
R=8.314;        %J/(mol.K) - gas constant
rho_cat=850;    %kg_cat/m3_reactor - density of the catalyst
Cpe=65;         %J/(mol.K) - Cp of ethylene
Me=28e-3;       %kg/mol - molar mass of ethylene
Meo=44e-3;      %kg/mol - molar mass of ethylene oxide
Hp=-2.11e5;     %J/mol O2 - enthalpy of Rp 
Hx=-4.42e5;     %J/mol O2 - enthaply of Rx
lG=0.04;        %W/(m.K) - conductivity of the reaction mixture
lC=0.22;        %W/(m.K) - conductivity of the used catalyst 
visc=1.5e-5;    %Pa.s - viscosity of the reaction mixture
eps=0.45;       %- - porosity (gas fraction) of the packed bed

%---------------------------------------------------------

%-------------------derived variables----------------------
Tin=Tin+273.2;  % changing Tin from C to K
Tc=Tin;         %K - coolant temp equals the inlet temp
C=P/(R*Tc);     %mol/m3 - molar concentration of the whole mixture == constant
CO2_0=fo2*C;    %mol/m3 - O2 con at z=0
CE_0=fe*C;      %mol/m3 - E con at z=0
%----tube diamater
if dtube == 1
    dt = 24.3e-3; %m - inner tube diameter in m
   
end
if dtube == 2
    dt = 49.3e-3;
   
end
if dtube == 3
    dt = 73.7e-2; 
    
end
%-------------------------------------------------------------



%------------- calculation domain & initial conditions
Lspan=[0:0.01:L]; %axial points at which output is generated
init=[CO2_0 CE_0 Tin]'; %initial vector for solving the set of odes



%----------------heat transfer model----------------------------
%------------- -------------------------------------------------
%the overall heat transfer coefficient is determined at the inlet
%conditions
REp=C*Me*dp*u/visc; %- - particle Reynolds number

Nu = 3.5*REp^0.7*exp(-4.6*dp/dt);
UU = Nu*lG/dt


%---------------------------------------------------------------

%---------------calling the ode solver--------------------------
options=odeset('RelTol',1e-6,'AbsTol',1e-9,'NormControl','on');
[z,x]=ode45(@der_eo,Lspan,init,options);
%---------------------------------------------------------------

%------------output variables-----------------------------------
%-------Axial profiles
n=20/0.01+1;                        %number of points in axial out variables
cono2=(CO2_0-x(2:n,1))/CO2_0;       %conversion of 02 as function of the reactor length
cone=(CE_0-x(2:n,2))/CE_0;          %conversion of E as function of the reactor length
ceo=-0.4*cono2*CO2_0+1.2*cone*CE_0; %concentration EO as function of the reactor length
s=1./((cone*CE_0)./ceo);            %EO selectivity based on E as function of the reactor length
temp=x(2:n,3);                      %temp as function of the reactor length
Tmax=max(temp)-273.2


%-------------plots--------------------------------------------
figure(1)
subplot(1,3,1)
plot(z(2:n),temp-273.2,'r')
xlabel('reactor length, m');
ylabel('temperature, C');


subplot(1,3,2)
plot(z(2:n),cono2,'r')
xlabel('reactor length, m');
ylabel('oxygen conversion, -');


subplot(1,3,3)
plot(z(2:n),s,'r')
xlabel('reactor length, m');
ylabel('selectivity (Seo)e, -');

%---------------------------------------------------------------





%---------------------------------------------------------------


%------------------------defining the differential equations------
function der=der_eo(z,x);

%x(1) = CO2
%x(2)=  CE
%x(3) = T

der=zeros(3,1);
kp=35.2*exp(-7200/x(3)); %mole o2 per kg cat per sec;
kx=74.1e3*exp(-10800/x(3)); %mole o2 per kg cat per sec;
effp=1; %- - effectiveness factor Rp, 
effx=1; %- - effectiveness factor Rx, 


der(1)=-(effp*kp+effx*kx)*x(1)*rho_cat/u; %mol balance O2: u*dCo2/dz=(-effp*kp*co2-effx*kx*co2)*rho_cat
der(2)=-(2*effp*kp+(1/3)*effx*kx)*x(1)*rho_cat/u; %mol balance E: u*dCE/dz=(-2*effp*kp*co2-1/3*effx*kx*co2)*rho_cat
der(3)=(-kp*effp*x(1)*rho_cat*Hp-kx*effx*x(1)*rho_cat*Hx-4*UU/dt*(x(3)-Tc))/(u*C*Cpe);
%energy balance:
%(u*C*Cp)*dT/dz=(-effp*kp*co2*Hp-effx*kx*co2*Hx)*rho_cat-4*U/dt*(T-Tc)
end

end
