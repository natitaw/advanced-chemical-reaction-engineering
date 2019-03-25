  function [ ]=FTsimple2016
% Simple model of a bubble slurry column for FT
% - isothermal
% - solubility co and h2 equal
% - rate: Rco = k*Cco(liquid)
% - large bubbles in plug flow
% - liquid - axially dispersed
% - small bubbles - axialy dispersed
% - only co balance is modelled
%
clear all
clc


% input parameters========================================================= 
P = 3e6;            % pa, pressure
T = 513;            % K, temp
Ug = 0.4;           % m/s overall superficial velocity gas
Ul = 0.01;          % m/s, overall superficial velocity liquid
L = 30;             % m reactor length
dt = 8;             % m reactor diameter
dct = 0.06;          % m colling tube diameter
dp = 50e-6;         % m diameter catalyst particle
h2co = 2;           % mole ration h2 over co inlet
ep = 0.25;          % volume fraction particles (catalyst) in the dense phase, m3_solids/m3_dense
Dsb = 10;           % dispersion coeff. small bubbles m2/s
Dl = 1;             % dispersion coeff. liquid
mh2 = 2.5;          % solubility parameter; m = Cg/Cl_interface
mco = 2.5;          % solubility parameter; m = Cg/Cl_interface
kr = 1.50e-4;       % reaction rate constant, mol Co / kg cat /s = overall rate, k = m3 / kgcat / s
rhocat = 2000;      % kg/m3_catalyst, density catalyst   
MWgas = 10e-3;      % average mol mass gas, kg/mol
rho_g_ref = 1.3;    % parameter in 
DHr = -152e3;        % J / mol / Co
HTC = 1000;          % W/m2/K 
DT = 15;            % temp difference between reactor and cooling

%==========================================================================

yco_in = 1/(1+h2co);    % mole craction co inlet
yh2_in = 1- yco_in;     % mole fraction h2 inlet
C =  P/(8.314*T)        % conc. mol/m3  
rho_g = MWgas*C ;       % density gas, kg/m3
Cco_in = yco_in*C       % conc co mol/m3
Ch2_in = yh2_in*C       % conc co mol/m3
A = 1/4*pi*dt^2;        % reactor diamter, m2
V = A*L;                % reactor violume, m3
Qv = A*Ug;              % volumetric flow rate inlet, m3/s
Act = pi*dct*L          % area for heat transfer of one cooling tube
Act1 = 1/4*pi*dct^2     % cylindrical area of one cooling tupe

CLH2 = Ch2_in/mh2;       % concentration H2 in liquid at equilibrium with inlet
CLCO = Cco_in/mco;      % idem
 

% hydrodynamics============================================================
V_small = 0.095*(1+0.8*ep/0.095)                    % rise velocity samll bubbles, m/s
e_DF = 0.6072*(1-0.7*ep/0.27)                       % hold-up gas in dense phase m3 gas / m3 dense phase
%ff = 0.27*(rho_g/rho_g_ref)^0.48
U_SB = V_small*e_DF                                 % superficial velocity small bubbles, m/s
U_LB = Ug - U_SB                                    % superficial velocity large bubbles, m/s 
e_LB = 0.7*U_LB^0.58                                % hold-up large bubbles, m3 gas / m3 reactor
e_gas = e_LB+e_DF*(1-e_LB)                          % total gas hold-up                     
%==========================================================================

% kla======================================================================
kla_LB = e_LB*2.25                                  % kla large bubble, 1/s, 4.5 is root of Dco/Dref, D is diffusion coefficient
kla_SB = e_DF*4.5                                   % kla small bubble, 1/s
%=========================================================================

% ========================================================================
% calling the solver
 solinit = bvpinit(linspace(0,L),[Cco_in Cco_in/mco  0 Cco_in 0]);
 sol = bvp4c(@f,@mat4bc,solinit);
 x = linspace(0,L)';
 y = deval(sol,x)';
 
 n=length(x);
 CLBout = y(n,1) % outlet concentration CO large bubbles
 CSBout = y(n,4) % outlet concentration CO small bubbles
 CLout = y(n,2)  % outlet concentration CO in the liquid.
 
  
 % calculation of the conversion
 FCOin = Ug*Cco_in
 
 FCOout = U_LB*CLBout+U_SB*CSBout+Ul*CLout
 
 XCo = 1- FCOout/FCOin
 
 % calculating the productivity
 P=FCOin*XCo*A/V                    %mol CH2/m3/s
 
 PP = P * 14e-3*V*3600*24/1000      %ton CH2/day   
 
 % heat removal 
 
 HR = P * V * abs(DHr)              % Watt
 
 WpCT = Act*HTC*DT                  % Watt per cooling tube
 
 NoCT = HR/WpCT                     % number of cooling tubes
 
 AreaRCT = NoCT*Act1/A
 
 figure(1)
 plot(x,y(:,1))
 xlabel('reactor length, m')
 ylabel('Co conc. in large bubbles, mol/m3')
 axis([0 L 0 300])
figure(2)
 plot(x,y(:,2))
 xlabel('reactor length, m')
 ylabel('Co conc. in liquid, mol/m3')
 axis([0 L 0 120])
 figure(3)
 plot(x,y(:,4))
  xlabel('reactor length, m')
 ylabel('Co conc. in small bubbles, mol/m3')
 axis([0 L 0 300]) 
 figure(4)
 plot(x,y(:,1)./y(:,2))
  xlabel('reactor length, m')
 ylabel('C large bubble / C Liquid')
 axis([0 L 0 5])
% ------------------------------------------------------------


% ------------------------------------------------------------



% CLB = y(1)
% CL = y(2)
% dCL/dx = y(3)
% CSB = y(4)
% dCSB/dx = y(5)


 
  function dydx = f(x,y)
    
     
        
     dydx = [-kla_LB*(y(1)/mco-y(2))/U_LB        
             y(3)
              (Ul*y(3) - kla_LB*(y(1)/mco-y(2))- kla_SB*(y(4)/mco-y(2))+(1-e_LB)*ep*rhocat*kr*y(2))/Dl
             y(5)
              (U_SB*y(5) + kla_SB*(y(4)/mco-y(2)))/Dsb]; 
  end
  
  % ------------------------------------------------------------
function res = mat4bc(ya,yb)
   
res = [  ya(1) - Cco_in
         ya(3) - Ul/Dl*(yb(2)-Cco_in/mco)
         yb(3)
         ya(5) - U_SB/Dsb*(ya(4)-Cco_in)
         yb(5)
         ];
 end
  




end  