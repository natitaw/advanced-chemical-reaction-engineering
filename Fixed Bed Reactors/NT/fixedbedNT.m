function[]=fixedbedNT

%1D fixed bed NT Shenkute

%input data
%what is the required data?
u=1;%superficial velocity m/s
a=10; %specific area
km=1;%mass transfer coefficient m/s
krxn=0.1;%reaction rate 1/s
n=2;%rxn order
Cf0=5; %inlet concentration mol/m3
C=40; %overall concentration mol/m3
L=1; %reactor length
Tin=20;%inlet temperature
Cp=40;
ap=200; %heat transfer coefficient W/m2-cat/K
DH=-1e4; %J/mol Enthalpy of reaction

%conditions at L=0; Start conditions

Cs0=km/(km+krxn)*Cf0; %To find the inlet concentration of catalyst

init=[Cf0
    Cs0
    Tin
    Tin];

%Empty matrix to solve eqn% This is how it works I guess

M=sparse(4,4); %goes in the odeset
M(1,1)=1;
M(3,3)=1; %Spots filled according to Matrix syntax

xspan=[0 L]; % span withwhich the differential will be solved dx


%Time to solve the ODE

options = odeset('Mass', M, 'MassSingular', 'yes', 'RelTol', 1e-6, 'AbsTol', 1e-6);

[z,c] = ode15s(@f, xspan, init', options); 

figure(1)
plot(z, c(:,1),'r', z, c(:,2))
xlabel('length, m');
ylabel('Conc. red= fluid, blue is solid');
figure(2)
plot(z,c(:,3),'r', z,c(:,4))
xlabel( 'length, m');
ylabel('T. red = fluid, blue is solid');

function out = f(z,c)

out = zeros(4,1);
out(1) =  -km*a/u*(c(1)-c(2));
out(2) = km*(c(1)-c(2))-krxn*c(2)^n;
out(3) = a*ap/(u*C*Cp)*(c(4)-c(3));
out(4) = ap*(c(4)-c(3))-krxn*c(2)^n*(-1*DH);

end
end


