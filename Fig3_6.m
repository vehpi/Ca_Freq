%Calcium Oscillation Frequency-Sensitive Gene Regulation and Homeostatic Compensation in 
%Pancreatic Beta-Cells. Bulletin of Mathematical Biology. 2017. 79-1295
%DOI: 10.1007/s11538-017-0286-1


% Matlab code for Oscillation Efficiency (Fig. 3) for 3
% different cooperativity values (nnE)
% Cooperativity vector

nnE=[1,2,4];

for i=1:3;
nE=nnE(i);

%Enzyme parameters
pE=0.09; dE=0.003;
D=10;
98
T=D:1:180;
KE=0:0.001:0.2;
[T,KE]=meshgrid(T,KE);

%Response to Oscillatory Ca
%Ebar represents the mean fractions of an activeted Enzyme
%During on times Ca is 0.1 and inbetween it is set to 0;
c0=0.1; %activation rate during on times

pEstar=pE*(c0^nE)./(c0^nE+KE.^nE);

%steady state fractions of active enzymes
E_ss=pEstar./(pEstar+dE);

%Activated enzyme franction (Eq. 12)
Ebar=E_ss.*(D./T+((pEstar./((pEstar+dE).*dE.*T)).*...
(1-exp(-D.*(pEstar+dE)))).*(1-exp(-dE.*(T-D)))...
./(1-exp(-(pEstar.*D+dE.*T))));

%Non-Oscillatory Ca with the same concentration
%subsrcipt c (_c) is used for indicating constant Ca
cc=c0*D./T;
E_c_ss=pE.*((cc.^nE)./(cc.^nE+KE.^nE))./(pE.*((cc.^nE)...
./(cc.^nE+KE.^nE))+dE);

%Oscillation Efficiency (Eq. 15)
Os_Ef=(Ebar-E_c_ss)./E_c_ss;

%Figure 3
figure(3)
subplot(3,1,i)
mesh(T,KE,Os_Ef)
end

% Matlab code for Fig. 6B,C

%Enzyme Parameters
nA=4; KcA=0.2; nI=4; KcI=0.2;
pA=0.09; dA=0.003; pI=0.09; dI=0.003;

%TF parameters
alphaA=0.01; betaI=0.01;

D=10;
T=D:1:180;
KA=0:0.001:0.28;
KI=0.4;
[T,KA]=meshgrid(T,KA);

%Response to Oscillatory Ca
%Abar and Ibar represent the mean fractions of active activator
%and inhibitor concentrations after sufficiently large
%oscillation cycles, respectively.

% During on times Ca is 0.1 and inbetween it is set to 0;
c0=0.1; %activation rates during on times

pAstar=pA*(c0^nA)/(c0^nA+KcA^nA);
pIstar=pI*(c0^nI)/(c0^nI+KcI^nI);

%steady state fractions of active enzymes
A_ss=pAstar./(pAstar+dA);
I_ss=pIstar./(pIstar+dI);

%Activated enzyme franctions
Abar=A_ss.*(D./T+((pAstar./((pAstar+dA).*dA.*T)).*...
(1-exp(-D.*(pAstar+dA)))).*(1-exp(-dA.*(T-D)))./...
(1-exp(-(pAstar.*D+dA.*T))));
Ibar=I_ss.*(D./T+((pIstar./((pIstar+dI).*dI.*T)).*...
(1-exp(-D.*(pIstar+dI)))).*(1-exp(-dI.*(T-D)))./...
(1-exp(-(pIstar.*D+dI.*T))));
Ainf=alphaA*(Abar.^2)./(KA.^2+Abar.^2);
Iinf=betaI*(Ibar.^1)./(KI.^1+Ibar.^1);
TFbar=Ainf./(Ainf+Iinf);

%Rates of Changs of infinity functions versus T
d_Abar_dT=A_ss*(-D./(T.^2)-(pAstar./(dA*(T.^2).*(pAstar+dA))).*...
((1-exp(-D.*(pAstar+dA)))).*(1-exp(-dA.*(T-D)))./...
(1-exp(-(pAstar.*D+dA.*T)))+...
(pAstar./((pAstar+dA).*dA.*T)).*...
(1-exp(-D*(pAstar+dA)))*dA.*(exp(-dA*(T-D))-...
exp(-(pAstar*D+dA*T)))./((1-exp(-(pAstar*D+dA*T))).^2));

d_Ibar_dT=I_ss*(-D./(T.^2)-(pIstar./(dI*(T.^2).*(pIstar+dI))).*...
((1-exp(-D.*(pIstar+dI)))).*(1-exp(-dI.*(T-D)))./...
(1-exp(-(pIstar.*D+dI.*T)))+...
(pIstar./((pIstar+dI).*dI.*T)).*...
(1-exp(-D*(pIstar+dI)))*dI.*(exp(-dI*(T-D))-...
exp(-(pIstar*D+dI*T)))./((1-exp(-(pIstar*D+dI*T))).^2));

d_Ainf_dT=2*(KA.^2).*(Abar.^...
(-2-1)).*d_Abar_dT./((1+(KA.^2).*(Abar.^-2)).^2);
d_Iinf_dT=1*(KI.^1).*(Ibar.^...
(-1-1)).*d_Ibar_dT./((1+(KI.^1).*(Ibar.^-1)).^2);
Diff=(d_Ainf_dT)-(d_Iinf_dT);
Sign_Diff=0.5*sign(Diff);

%Figure 6 panels B and C
figure()
subplot(1,2,1)
mesh(T,KA,TFbar)
subplot(1,2,2)
mesh(T,KA,Sign_Diff)