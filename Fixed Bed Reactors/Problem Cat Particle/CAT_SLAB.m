function []=CAT_SLAB() 
% function solves catalyst particle. slab geometry. first order reaction.
% Interal and external mass transfer.
% equation:
% c'' - lamda^2*c = 0; c'' = second derivative of c to theta
% lamda = L*(kv/De)^0.5 = Thiele modulus
% c is the dimensionless concentration c/c in bulk
% theta is dimensionless length of the slab x/L, 0 is the centre line, 1 is
% outside
% boundary conditions:
% at the centre line c'=0
% at the outside c'=Bi*(1-c), c' = first derivative of c to theta
% BVP

clc
clear all

% input parameters
lamda = 1;                                  % value of lamda.. value can be changed. 
Bi = 1;                                     % value of Bi.. value can be changed. 

n=100; % number of points on the theta axis. post-processed; the solver selects and adapts the grid size

function der = der_y(theta,y)               % function to solvbed by BVP4C

der=[y(2)                                   % y(1) = C; y(2) = dy(1)/dtheta
    lamda^2*y(1)];                          % dy(2)/dtheta = lamda^2 * y(1)....main Differential Equation   
end

function res=bc(y0,y1)                      % Boundary conditions.. y0 is left boundary, i.e theta = 0 ; y1 is right boundary, i.e theta = 1
res=[y0(2)                                  % Index (1) or (2) means y(1) which is  c or y(2)which is dcdtheta
    y1(2)-(Bi*(1-y1(1)))];                  % Boundary conditions equations..BC1...y0(2) = 0.
end                                         % BC2 y1(2) - Bi(1-y1(1)) = 0



theta = linspace(0,1,n+1)';                 % defining the theta grid

sol=bvpinit(theta,[1 0.5]);                 % Initial guess for  [y(1) y(2)]; here a constant value for y(1) and y(2) is used (same iniitial valua at all grid points.  
                                            
sol=bvp4c(@der_y,@bc,sol);                   % Boundary value solver, bvp4c (derivative, boundary condition, solution points)
C = sol.y';                                  % C has 2 columns... column 1 is Concentration and column 2 is dConcentration/dtheta

  
Ca = cosh(lamda.*theta)./(cosh(lamda)+lamda/Bi*sinh(lamda)); % elementary solution 

% effectiveness factor

eff = tanh(lamda)/lamda/(1+lamda/Bi*tanh(lamda)) % elementary effectiveness factor

effn= trapz(theta,C(:,1))                        % numerical effectivenes factor


plot (theta,C(:,1),'ro',theta,Ca)                % plot of Concentration against theta
xlabel('dimensionless length, -');
ylabel('dimensionless Conc. of A, -');
end