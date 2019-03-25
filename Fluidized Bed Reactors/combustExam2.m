function [ ] = combustExam2(Ht,u0,L,dt,T,Cb) 
% code to dimension a FCC combustor 
% algorithm:
% -solve O2 balance with gues value for Cc (mol C per kg / catalyst)
% -compute the coke balance
% -if the coke balance is satisfied keep the solution if not
% -new Cc value until coke balance is satisfied
% for the iteration Fminsearch is used
clc
clear all

% input paramters
% ==============================================
u0 = 0.95  ;
L = 5 ;
dt = 14;


T = 900;
%==============================================



Ht = L^0.5*dt^0.42;

Cb = 2.8; % concentration of O2 at the inlet


rho_cat = 500; % kg/m3 reactor --> estimated from 1500 kg/m3 particle * porosity of 2/3
k=1e9*exp(-18889/T);


% performing the iteration":

Ccinit = 0.45; % guess value of Cc

solfc = fminsearch(@balance,Ccinit);

Cc = solfc;
Xc=1-Cc                                                 % coke conversion

% calculation the solution with the correct Cc value

[x,O2b,O2d,balC]=FBDeemterComExam2(Ht,u0,L,dt,T,Cb,Cc);

n=length(x);
oxygen_fed = 1/4*pi*dt^2*u0*1e5/(8.314*T)*0.21*32e-3    % oxygen fed in mol/s
Ofedm = oxygen_fed/0.23*3600/1000                       % oxygen in ton/h
CO2d = mean(O2d)                                        % average O2 concentration in the dense phase
coke_in = 600*1                                         % inflow of coke
coke_out = 600*Cc                                       % outflow of coke
coke_reacted = (4/3)*Cc*rho_cat*k*CO2d*L*1/4*pi*dt^2    % coke reacted in CSTR
oxygen_reacted = (Cb-O2b(n))*1/4*pi*dt^2*u0
oxygen_conversion = 1-O2b(n)/Cb                         % oxygen conversion

balC=coke_in-coke_out-coke_reacted

plot(x,O2b,x,O2d)




% file defining the coke balance
    function f = balance(Cc)


        
        [x,O2b,O2d,balC]=FBDeemterComExam2(Ht,u0,L,dt,T,Cb,Cc);


        f = abs(balC);

    end



end

