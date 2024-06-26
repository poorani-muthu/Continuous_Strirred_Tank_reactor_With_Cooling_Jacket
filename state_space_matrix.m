
syms Ca T Tc q
A=jacobian(plantjac1,[Ca,T]);
B=jacobian(plantjac1,[Tc,q]);
C=[1 0;0 1];
D=0;
Ca=0.989;
T=296.616;
Tc=270;
q=100;
AA=eval(A);
BB=eval(B);
reac_cont=ss(AA,BB,C,D);
reac_discjac=c2d(reac_cont,0.1);

function dCdt = plantjac1()
syms T Ca Tc q
V=100;
EoverR=8750;
rho=1000;
Cp=0.239;
mdelH=5*(10^4);
ko=7.2*(10^10);
Caf=1;
Tf=350;
UA=50000;
dCdt(1)=(q*(Caf-Ca)*(1/V))-(ko*exp(-EoverR/(1*T))*Ca);
dCdt(2)=((q*rho*Cp*(Tf-T))+(mdelH*V*ko*exp(-EoverR/(1*T))*Ca)+(UA*(Tc-T)))/(V*rho*Cp);
end