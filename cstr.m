                
% Open Loop Simulation of CSTR
Tc=270;
q=100;

tspan=[0 5];
cO=[0.9, 305];

[t,c]=ode15s(@(t,C)cstr_1(t,C,q,Tc),tspan,cO);

figure(1)
sgtitle("Openloop Simulation of CSTR with JC")
subplot(2,1,1)
plot(t,c(:,2),'--b')
title('Temperature Vs Time')
subplot(2,1,2)
plot(t,c(:,1),'--r')
title("Reactant concentration Vs Time")


% CSTR FUNCTION FILE %

function dCdt = cstr_1(t,C,q,Tc)

    dCdt=zeros(2,1);

    %Constants%

    V=100; %m3
    EoverR=8750;
    rho=1000; %kg/m3
    Cp=0.239;%J/kg K
    mdelH=5*(10^4);%J/mol
    ko=7.2*(10^10); %1/sec
    Caf=1; %m3/mol
    Tf=350;%K
    UA=50000;

    %Equations%

    dCdt(1)=(q*(Caf-C(1))*(1/V))-(ko*exp(-EoverR/(1*C(2)))*C(1));
    dCdt(2)=((q*rho*Cp*(Tf-C(2)))+(mdelH*V*ko*exp(-EoverR/(1*C(2)))*C(1))+(UA*(Tc-C(2))))/(V*rho*Cp);
    
end