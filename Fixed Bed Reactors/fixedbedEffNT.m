function fixedbedEffNT
clc
clear all

u=1;
km=0.001;
kv=0.2;%reaction rate
n=1.5; %reaction order
Cf0=5; %inlet Cf
L=10; %Reactor length
Lp=2e-3; %particle length
eps = 0.5; %porosity
a=6/Lp*(1-eps);; %specific area
aa = 6/Lp; %specific area m2 external cat surface/m3 particle
D=1e-5;%m2/s

%initial conditions at L=0
 Cs0=Cf0;
 

init = [Cf0
    Cs0];

M = sparse(2,2);
M(1,1)=1;

xspan=[0:0.01:L];

options = odeset ('Mass', M, 'MassSingular', 'yes', 'RelTol', 1e-6, 'AbsTol', 1e-6);

[z c] = ode15s(@f, xspan, init', options);

p=length(z);

for i = 1:p
    effp(i)=CAT_SLAB_nth(kv, D, Lp, c(i,2),n);
    efft(i)=aa*km*(c(i,1)-c(i,2))./(kv*c(i,1).^n); %I don't get this part??
end
figure(1)
subplot(1,2,1)
plot(z,c(:,1),'r', z, c(:,2),'k')
xlabel(' length, m');
ylabel('Conc. red = fluid, black is solid');
subplot(1,2,2)
plot(z,effp,z,efft,'r')
xlabel(' length, m');
ylabel('eff. factor. red = total, blue is intera particle');
function out = f(z,c)

out = zeros (2,1);

[eff]=CAT_SLAB_nth(kv,D,Lp,c(2),n) %calling the function
out(1)=-km*a/u*(c(1)-c(2));
out(2)=-kv*c(2)^n*eff+km*aa*(c(1)-c(2));
end 
end 
