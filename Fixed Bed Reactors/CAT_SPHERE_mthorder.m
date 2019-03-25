function []=CAT_SPHERE_mthorder() 
% function solves catalyst particle. slab geometry. first order reaction.
% Interal and external mass transfer.
% equation:
% c'' + 2/theta*c' - lamda^2*c^m = 0; cii = second derivative of c to theta.
% lamda = R*(kv/De)^0.5
% c is the concentraton 
% theta is dimensionless radius of the sphere r/R, 0 is the centre line, 1 is
% outside
% boundary conditions:
% at the centre line c'=0
% at the outside c'=Bi*(Cb-c)
% Cb is the concentration in the bulk


clc
clear all


% setting the parameters
R = 1e-3;                                   % radius of the paricle
De = 1e-6;                                  % eff. diffusion coef
Kv = 5;                                   % reaction rate constant Rate mol/m3 particle /s
lamda = R*(Kv/De)^0.5                      % value of lamda.. value can be changed. Same for Bi
Bi = 100000
m = 2;                                      %  reaction order       
Cb = 10;                                    %  bulk concentration



n=100; % number of points on the theta axis. post-processed; the solver selects and adapts the grid size

function der = der_y(theta,y)              % derivative of y wrt theta

der=[y(2)                                  % y(1) = C; y(2) = dy(1)/dtheta
    lamda^2*y(1)^m];                       % dy(2)/dtheta = lamda^2 * y(1)....main Differential Equation , the term 2/thetay(2) we put in via the S matrix  
end

function res=bc(y0,y1)                      % Boundary conditions.. y0 is left boundary, i.e theta = 0 ; y1 is right boundary, i.e theta = 1
res=[y0(2)                                  % Index (1) or (2) means y(1) which is c or y(2)which is dcdtheta
    y1(2)-(Bi*(Cb-y1(1)))];                 % Boundary conditions equations..BC1...y0(2) = 0...
end                                         % BC2.  y1(2) - Bi(Cb-y(1)) = 0



theta = linspace(0,1,n+1)';                 % defining the theta grid

solinit=bvpinit(theta,[0 0.5]);             % Initial guess for the bvp4c. theta axis... min = 0, max = 1, divisions = 101.. the more the divisions the smoother the curve..

S = [0 0 ; 0 -2];                           % defining S                 

options =bvpset('SingularTerm',S);



sol=bvp4c(@der_y,@bc,solinit,options);                  
C = sol.y';      % C has 2 columns... column 1 is Concentration and column 2 is dConcentration/dtheta

%calculating the effectiveness factor 
j = length(C);

eff = 3*De*C(j,2)/(Kv*Cb^m*R*R)     % effectiveness factor =flux divided by reaction rate at bulk conditions. 3/R comes from area over volume ratio, 
%extra R is to undo the dimenensionelessness of dCdtheta

% analytical solution for effectiveness factor for thiele > 2.5 and Bi >> to
% check our numerical work

Thieleb = lamda/3*((m+1)/2*Cb^(m-1))^0.5
effb = 1/(3*Thieleb^2)*((3*Thieleb)*coth(3*Thieleb)-1)


plot (theta,C(:,1),'ro')           % plot of Concentration against theta
xlabel('dimensionless length, -');
ylabel('Conc. of A, -');
end