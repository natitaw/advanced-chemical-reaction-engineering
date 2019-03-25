function [x,O2b,O2d,balC]=FBDeemterComExam2(Ht,u0,L,dt,T,Cb,Cc)
%FIRST FILE for question with vandemeter model
% function solves the oxygen mass balance using the van deemter model
% Cb = C in bubble, Cd = C in dense
% derivaties to theta (x/L) demennsioneless lenght of the reactor
% Cb' = - Nt*(Cb-Cd), Cb' first derivative to theta
% Cd'' =  -Nt*Ne*(Cb-Cd) + Nr*Nr*Cd, Cd'' is second derivative to theta
% SRA Kersten 2015                                    


% parameters:
D=1e5; % dispersion coefficient; very large to ensure CSTR of dense phase

rho_cat = 500; % kg/m3 reactor --> estimated from 1500 kg/m3 particle * porosity of 2/3

k=1e9*exp(-18889/T); % given reaction rate 

Cbin=Cb; % concentration of Oxyegn at the inlet (T = 900 K)

% Nt, Ne, Nr calculation
% Cc; mol C per kg catalyst is quess value.
Nt  = L/Ht;
Nr  = Cc*rho_cat*k*L/u0;
NE  = u0*L/D;


%------------------ solving the equations -------------------
solinit = bvpinit(linspace(0,1),[1 0 1]);
sol = bvp4c(@f,@mat4bc,solinit);
x = linspace(0,1)';
y = deval(sol,x)';

% checking the solution: 
grad = y(:,3)-y(:,1);
a=trapz(x,Nr*y(:,1)); % reaction in dense phase
aa=trapz(x,Nt*grad);  % transfer from bubble to dense
check = aa-a;         % transfer = reaction

% computing the coke balance
O2b=y(:,3);
O2d=y(:,1);
CO2d = mean(O2d);
coke_in = 600*1;
coke_out = 600*Cc;
coke_reacted = (4/3)*Cc*rho_cat*k*CO2d*L*1/4*pi*dt^2;

balC=coke_in-coke_out-coke_reacted;

%plot(x,O2b,'r',x,O2d)



% ------------------------------------------------------------


% ------------------------------------------------------------
% funtions (ODE and the boundary conditions
% Cdense = y(1)
% dCdense/dx = y(2)
% Cbubble = y(3)
 
  function dydx = f(x,y)
     dydx = [            y(2) 
             -Nt*NE*(y(3)-y(1))+Nr*NE*y(1)
             -Nt*(y(3)-y(1))]; 
  end
  
% ------------------------------------------------------------
function res = mat4bc(ya,yb)
res = [  ya(2)
         yb(2) 
         ya(3)-Cbin];
end
  

end  