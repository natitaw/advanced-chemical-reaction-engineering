function []=fixedbedNT2

%1D fixed bed
u=1;
a=10;
km=1;
krxn=0.1; %reation rate
n=2;%rxn order
rhof=1; %density of f
DH=-1e4; %enthalpy of rxn
Cf0=5; %inlet concentration of f
C=40; %toatal concentration
L=1; %reactor length
Tin=20;%inlet T
Cp=40;
ap=200;%alpha p


Cs0=km/(km+krxn)*Cf0;

init=[Cf0
    Cs0
    Tin
    Tin];

M=sparse(4,4);
M(1,1)=1;
M(3,3)=1;

xspan=[0 L];

options = odeset('Mass', M, 'MassSingular', 'yes', 'RelTol', 1e-6, 'AbsTol', 1e-6);

[z,c]=ode15s(@f, xspan, init', options);
figure(1)
plot(z,c(:,1),'r', z, c(:,2))
xlabel('length, m');
ylabel('Conc. red = fluid, blue is solid');
figure(2)
plot(z,c(:,3),'r', z,c(:,4))
xlabel('length, m');
ylabel('T. red=fluid, blue is solid');

function out = f(z,c)

out = zeros(4,1);
out(1)=-km*a/u*(c(1)-c(2));
out(2)=km*(c(1)-c(2))-krxn*(c(2))^n;
out(3)=a*ap/(C*Cp*u)*(c(4)-c(3));
out(4)=ap*(c(4)-c(3))-krxn*(c(2))^n*(-1*DH);
end
end

