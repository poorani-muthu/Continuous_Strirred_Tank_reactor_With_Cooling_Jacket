clear all
% PID controller for CSTR with cooling jacket
% manipulated variables : q and Tc
% controll variables : Ca and T
% Ca is controlled by q
% T is controlled by Tc

Kpc= 0.1;
tauic=50;
taudc=2.5;

KpT =0.1;
tauiT =50; 
taudT = 2.5;

time=0:0.1:30;
delt=0.1;

init(1,:)=[0.98,290]; %initial concentraion and temerarutre of outlet stream
spT=296.615795335250; %set point for temp
spCa=0.988382968500109; %setpoint for height

%initial error
eT=abs(init(1,2)-spT);
eCa=abs(init(1,1)-spCa);

Cmatf(1,:)=init(1,:);
q_(1)=100; 
Tc_(1)=270;
TF(1)=0;

for i=2:length(time)
    q_new=q_(i-1);
    Tc_new=Tc_(i-1);

    t(i,:)=[time(i-1) time(i)]; %loop runs for every 0.1 sec and last value is stored in Cmatf

    [tmat,Cmat]=ode45(@(t,C)CSTR_J(t,C,q_new,Tc_new),t(i,:),init);

    TF(i)=tmat(end);

    Cmatf(i,:)=Cmat(end,:); 
    init=Cmat(end,:); %initial cond changed to previous iteration values

    eT(i)=abs(Cmatf(i,2)-spT); %error b/w sp and cstr temp at different time
    eCa(i)=abs(Cmatf(i,1)-spCa);

    q_(i)=q_(1)+(Kpc*(eCa(i)+((delt/tauic)*(sum(eCa)))+((taudc/delt)*(eCa(i)-eCa(i-1)))));
    Tc_(i)=Tc_(1)+(KpT*(eT(i)+((delt/tauiT)*(sum(eT)))+((taudT/delt)*(eT(i)-eT(i-1)))));
    vec = [i Tc_(i) TF(i)]
end

Canew=Cmatf(:,1);
Tnew=Cmatf(:,2);

tplot=0:0.1:TF(end);
Casp=spCa.*ones(1,length(time));
Tsp=spT.*ones(1,length(time));

%plotting
figure(1)
sgtitle("PID Control on MIMO Model of CSTR")

subplot(2,2,1)
plot(tplot,Casp,'k')
hold on
plot(tplot,Cmatf(:,1),'b')
xlabel('time sec')
ylabel('Ca')
title('Concentration and Set point')
ylim([0.95 1])

subplot(2,2,2)
plot(tplot,Tsp,'k')
hold on
plot(tplot,Cmatf(:,2),'b')
xlabel('time sec')
ylabel('T')
title('Temperature and Set point')
ylim([290 300])

subplot(2,2,3)
plot(tplot,q_,'r')
xlabel('time sec')
ylabel('m^3/s')
title('Manipulated variable wh Vs. time')
ylim([97 103])

subplot(2,2,4)
plot(tplot,Tc_,'r')
xlabel('time sec')
ylabel('Tc')
title('Manipulated variable wc Vs. time')
ylim([260 275])

hold off

function dCdt = CSTR_J(t,C,q,Tc)
    
    %q=100; %m3/sec
    V=100; %m3
    EoverR=8750;
    rho=1000; %kg/m3
    Cp=0.239;%J/kg K
    mdelH=5*(10^4);%J/mol
    ko=7.2*(10^10); %1/sec
    Caf=1; %m3/mol
    Tf=350;%K
    UA=50000;

    dCdt(1,1)=(q*(Caf-C(1))*(1/V))-(ko*exp(-EoverR/(1*C(2)))*C(1));
    dCdt(2,1)=((q*rho*Cp*(Tf-C(2)))+(mdelH*V*ko*exp(-EoverR/(1*C(2)))*C(1))+(UA*(Tc-C(2))))/(V*rho*Cp);
end