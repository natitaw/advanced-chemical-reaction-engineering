 function [eff]=CAT_SLAB_nth(kv,D,L,Cs,m) 
% Function solves catalyst particle. slab geometry. nth  order reaction.
% Only Interal mass transfer.
% equation:
% c'' - lamda^2*c = 0; cii = second derivative of c to theta.
% lamda = L*(kv/De)^0.5 
% c is the concentration
% theta is dimensionless length of the slab x/L, 
% 0 is the centre line, 1 is outside
% boundary conditions:
% at the centre line c'=0
% at the outside 
% the effectiveness factor is calculated
% input: 
% kv, reaction rate on volumetric basis, mol / m3 catalyst / s
% D, effective difusion coefficient, m2/s
% L, characteristic length, m
% Cs, concentration at the surface, mol/m3
% m order of the reaction



kv = 10;
D = 1e-5;
L=2e-3;
Cs = 1;
m = 2; %reaction order

lamda = L*(kv/D)^0.5;




n=100; % number of points on the theta axis. post-processed; the solver selects and adapts the grid size

function der = der_y(theta,y)               % derivative of y wrt theta

der=[y(2)                                   % y(1) = C; y(2) = dy(1)/dtheta
    lamda^2*y(1)^m];                        % dy(2)/dtheta = lamda^2 * y(1)....main Differential Equation   
end

function res=bc(y0,y1)                      % Boundary conditions.. y0 is left boundary, i.e theta = 0 ; y1 is right boundary, i.e theta = 1
res=[y0(2)                                  % Index (1) or (2) means y(1) which  c or y(2)which is dcdtheta
    y1(1)-Cs];                              % Boundary conditions equations..BC1...y0(2) = 0....@ yo, y(2) = 0 
end                                         % BC2..... y1(2) - (Bi*(1-y1(1)))..... @ y1, y(2) - Bi(1-y(1)) = 0



theta = linspace(0,1,n+1)';

sol=bvpinit(theta,[1 0.5]);                 % Grid,  Initial guess for  [y(1) y(2)]; here a constant value for y(1) and y(2) is used (same iniitial valua at all grid points.  
                                            
sol=bvp4c(@der_y,@bc,sol);                  % Boundary value solver, bvp4c (derivative, boundary condition, solution points)
C = sol.y';                                 % C has 2 columns... column 1 is Concentration and column 2 is dConcentration/dtheta

mm=size(C);

flux = D*C(mm(1),2)/L;

rbulk = L*kv*Cs^m;

eff = flux/rbulk;


end