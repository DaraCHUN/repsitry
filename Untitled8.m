clc;
clear;
f = linspace(0.3,5,1000); 
c=3e8;
lambda=c./f.*1e-9;
G1= 4.597653;
G2=24.683648;
d1=0.1224987 ;
d2=0.0730035;
Z0 = 377;
%
ep1=(1+0.039.*G1.*f.^(-0.4)-j.*0.03.*G1.*f.^(-0.4)).^2;
ep2=(1+0.039.*G2.*f.^(-0.4)-j.*0.03.*G2.*f.^(-0.4)).^2;
gamma1=j.*2.*pi.*sqrt(ep1)./lambda;
gamma2=j.*2.*pi.*sqrt(ep2)./lambda;
%
t1= tanh(gamma1.*d1);
t2= tanh(gamma2.*d2);
%
Zc1=Z0./sqrt(ep1);
Zc2=Z0./sqrt(ep2);
%
Z2=Zc2.*t2;
Z1=Zc1.*((Z2+Zc1.*t1)./(Zc1+Z2.*t1));
%
Gamma=(Z0-Z1)./(Z0+Z1);
Gamma_abs = abs(Gamma);
Gamma_dB = 20*log10(Gamma_abs);
%
plot(f,Gamma_abs,'b','linewidth',2');
xlabel('f(GHz)');
ylabel('|\Gamma|');
grid on;