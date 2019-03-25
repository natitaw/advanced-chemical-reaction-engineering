function[yCOeq, yH2eq, yCO2eq, yH2Oeq, yCH3OHeq]=Equilibrium_model(yCH3OH, yCO, yCO2, yH2,yH2O, T, P)
R=8.314; %Gas constant in J K-1 mol
%-----------------------Solver--------------------------------
X0=[0.01 0.01]'; 
options=optimset('TolFun',1e-9,'TolX',1e-6, 'Display','off');
X = fsolve(@fun,X0,options);

    function F=fun(X)
%X(1) conversion CH3OH, X(2) conversion H2O 
 

Ka1=exp((1/(R.*T)).*(7.44140e4+1.89260e2.*T+3.2443e-2.*T.^2+7.0432e-6.*T.^3-5.6053e-9.*T.^4+1.0344e-12.*T.^5-6.4364e1.*T.*log(T)));
Ka2=exp((1/(R.*T)).*(-3.94121e4-5.41516e1.*T-5.5642e-2.*T.^2+2.5760e-5.*T.^3-7.6594e-9.*T.^4+1.0161e-12.*T.^5-1.8429e1.*T.*log(T)));

[fc_CH3OH, fc_CO, fc_CO2, fc_H2, fc_H2O]=SKR_model(T, P); 
 
%(1)CO+H2=CH3OH  
%(2)CO2+H2=CO+H2O

Denom=(1-2.*X(1).*yCO);
yCOeq= (yCO-X(1).*yCO+yCO2.*X(2))./Denom; 
yH2eq=(yH2-2.*yCO*X(1)-yCO2.*X(2))./Denom;
yCO2eq= (yCO2-yCO2.*X(2))./Denom;
yH2Oeq= (yH2O+X(2).*yCO2)./Denom;
yCH3OHeq= (yCH3OH+X(1).*yCO)./Denom;


Kp1= (yCH3OHeq./(yCOeq.*yH2eq.^2)).*(1./P.^2);
Kf1= (fc_CH3OH)./(fc_CO.*fc_H2.^2); %Methanol synthesis
Kp2= (yCOeq.*yH2Oeq)./(yCO2eq.*yH2eq);
Kf2= (fc_CO.*fc_H2O)./(fc_CO2.*fc_H2); %Water-gas shift reaction

F=zeros(2,1); 
F(1)= Ka1-Kp1.*Kf1;
F(2)= Ka2-Kp2.*Kf2;
    end 
end