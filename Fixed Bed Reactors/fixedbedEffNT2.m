function fixedbedEffNT2
clc
clear all
u=1;
eps=0.5;
km=0.001;
n=1.5;
kv=0.2;
Lp=2e-3;
a=6*(1-eps)/Lp;
Cf0=5;
L=10;

aa=6/Lp;
D=1e-5;


Cs0=Cf0;

init = [Cf0
    Cs0];

M=sparse(2,2);
M(1,1)=1;

xspan=[0:0.01:L];

options = odeset('Mass', M, 'MassSingular', 'yes', 'RelTol', 1e-6);

[z c]=ode15s(@f, xspan, init',options);

p=length(z);

for i=1:p
    effp(i)=CAT_SLAB_nth(kv, D, Lp, c(i,2),n);
    efft(i)=aa*km*(c(i,1)-c(i,2))./(kv*c(i,1).^n);
end
figure(1)
subplot(1,2,1)
plot(z,c(:,1),'r', z, c(:,2), 'k')
xlabel('length,m');
ylabel('Conc. red - fluid, black is solid');
subplot(1,2,2)
plot(z,effp,z,efft, 'r')
xlabel('length, m');
ylabel('eff. factor. red=total, blue is intera particle');
function out=f(z,c);

out = zeros(2,1);
[eff]=CAT_SLAB_nth(kv, D, Lp, c(2),n);
out(1)=-km*a/u*(c(1)-c(2));
out(2)=kv*(c(2))^n*eff-km*aa*(c(1)-c(2));
end
end
