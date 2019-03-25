function [X, Xpfr, Xcstr]=FBdeemterL_stripped
% This code solves the Van Deemter BFB model for a first order reaction and
% under the assumption that there is no flow going through the dense phase.
% The following euguations are solved:
% Cb' = - Nt*(Cb-Cd)
% Cd'' = -Nt*NE*(Cb-Cd) + Nr*NE*Cd
% Cb is the concentration of the reactant in the bubble phase
% Cd is the concentration of the reactant in the dense phase
% Cb' is the first derivative of Cb to the dimensionless height (x)
% Cd'' is the second derivative of Cd to the dimensionless height (x)
% Boundary conditions:
%   x = 0 --> Cb = Cb_in, Cd' = 0
%   x = 1 --> Cd' = 0
% The function calculates the conversion (X), and the conversions of the
% equivalent PFR (Xpfr) and CSTR (Xcstr). Equivalent means equal Nr.

% S.R.A. Kersten 2016


% definition of Nt, Nr, NE
Nt=2;
Nr=1;
NE=2;


% Cbin 
Cbin=1;


% calling the solution
solinit = bvpinit(linspace(0,1),[1 0 1]);
sol = bvp4c(@f,@mat4bc,solinit);
x = linspace(0,1)';
y = deval(sol,x)';

% plotting the solution
% Cdense = y(1)
% dCdense/dx = y(2)
% Cbubble = y(3)
 plot(x,y(:,1),'b',x,y(:,3),'r')
 ylabel('concentration, blue is dense phase, red is bubble phase')

% calculation of the outlet concentration (bubble phase)
n=length(y);

cout=y(n,3);

X = 1-cout                 % conversion
coutfixed = Cbin*exp(-Nr) ;
CoutCSTR = Cbin/(1+Nr);
Xcstr = 1- CoutCSTR        %  conversion CSTR
Xpfr = 1 - coutfixed       % conversion  pfr
eta = X/Xpfr               % gas-solid contacting
% definition of the function (dydx) & the boundary conditions (res)

 
    function dydx = f(x,y)
    % Cdense = y(1)
    % dCdense/dx = y(2)
    % Cbubble = y(3)
        
        
        dydx = [   y(2)     
            -Nt*NE*(y(3)-y(1))+Nr*NE*y(1)
               -Nt*(y(3)-y(1)) ]; 
    end
  

    function res = mat4bc(ya,yb)
    
        res = [ ya(2)
            yb(2)
              ya(3)-Cbin ];
    end
 

end  