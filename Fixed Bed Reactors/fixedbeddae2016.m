function  []=fixedbeddae2016
%1D fized bed reactor with surface reacion
% non-isothermal

%input data
%--------------------------------------------------------------------

u=1; % superficial velocity m/s
a=10; % specific area m2 external cat surface / m3 reactor
km=1; %  mass transfer coefficient m/s
kr=0.1; % reaction rate 1/s
n = 2; % order
inletCf=5; % inlet concentration mol/m3
C = 40; % overall concentration
L = 1; % reactor length
Tin = 20; % inlet temperature
Cp = 40; % J/mol/K
alpha = 200; % W/m2-cat/K 
DH = -1e4; % J/mol

%  conditions at L =0;
inletCs=km/(km+kr)*inletCf;

init=[inletCf
      inletCs
      Tin
      Tin];



% A constant, singular mass matrix
M = sparse(4,4);
M(1,1) = 1;
M(3,3) = 1; 

%

xspan = [0 L];

% Use the LSODI example tolerances.  The 'MassSingular' property is

options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-6,'AbsTol',1e-3);

[z,c] = ode15s(@f,xspan,init',options);

%analytical solution first order reaction:

%Ko = (1/km+1/kr)^-1;
%ca=inletCf*exp(-a*Ko.*z/u);

figure(1)
plot(z,c(:,1),'r',z,c(:,2))
xlabel(' length, m');
ylabel('Conc. red = fluid, blue is solid');
figure(2)
plot(z,c(:,3),'r',z,c(:,4))
xlabel(' length, m');
ylabel('T. red = fluid, blue is solid');
% --------------------------------------------------------------------------

function out = f(z,c)

% c(1)  = Cf
% c(2) = Cs
% c(3) = Tf
% c(4) = Ts

out=zeros(4,1);
out(1) = -km*a/u*(c(1)-c(2));
out(2) = km*(c(1)-c(2))-kr*c(2)^n;
out(3) = a*alpha/(u*C*Cp)*(c(4)-c(3));
out(4) = alpha*(c(4)-c(3))-kr*c(2)^n*(-1*DH);
end


end