function  []=fixedbed_with_effectiveness %doesn't have to be a function, can be a script
%1D fized bed reactor with surface reacion

clear all
%input data
%--------------------------------------------------------------------

u=1; % superficial velocity m/s
km=0.001; %  mass transfer coefficient m/s
kv=0.2; % reaction rate giving the rate in mol / m3 cat /s
m = 1.5; % order
inletCf=5; % inlet concentration mol/m3
L = 10; % reactor length, m
Lp = 2e-3;  % particle L, m
D = 1e-5; % m2/s
eps = 0.5; % porosisty, m3 fluid / m3 reactor
a=6/Lp*(1-eps); % specific area m2 external cat surface / m3 reactor
aa = 6/Lp ; % pecific area m2 external cat surface / m3 particle

%  conditions at L =0;
inletCs=inletCf;

init=[inletCf
      inletCs];
%BC at x=0


% A constant, singular mass matrix
M = sparse(2,2); %better way to save the 0's 
M(1,1) = 1;
 

%

xspan = [0:0.01:L];

% Use the LSODI example tolerances.  The 'MassSingular' property is

options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-6,'AbsTol',1e-6);

[z,c] = ode15s(@f,xspan,init',options);

p=length(z);

for i = 1:p
    effp(i)=CAT_SLAB_nth(kv,D,Lp,c(i,2),m);
    efft(i)=aa*km*(c(i,1)-c(i,2))./(kv*c(i,1).^m); %post process to get eff
end




%analytical solution first order reaction:

%Ko = (1/km+1/kr)^-1;
%ca=inletCf*exp(-a*Ko.*z/u);

figure(1)
subplot(1,2,1)
plot(z,c(:,1),'r',z,c(:,2),'k')
xlabel(' length, m');
ylabel('Conc. red = fluid, black is solid');
subplot(1,2,2)
plot(z,effp,z,efft,'r')
xlabel(' length, m');
ylabel('eff. factor. red = total, blue is intera particle');
% --------------------------------------------------------------------------

function out = f(z,c) %derivative

% c(1)  = Cf
% c(2) = Cs

[eff]=CAT_SLAB_nth(kv,D,Lp,c(2),m); %calling the function
%eff = 1

out=zeros(2,1);
out(1) = -km*a/u*(c(1)-c(2));
out(2) = aa*km*(c(1)-c(2))-eff*kv*c(2)^m;
%the two equations
end


end