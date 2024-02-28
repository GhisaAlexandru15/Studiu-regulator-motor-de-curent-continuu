Ra=0.92;
La=10^(-3);
Ke=0.296;
Kt=0.294;
J=7*10^(-4);
Bf=3.35*10^(-4);
KPWM=38.46;
A=[-Ra/La, -Ke/La, 0; 
    Kt/J , -Bf/J , 0; 
    0    , 1     , 0];
B=[KPWM/La ; 0 ; 0];
C=[ 0 0 1];
D=[0];
Co=ctrb(A,B)
rank(Co)
p1=-89.9132
p2=-53.9855
K=acker(A,B,[p1, p1,p1])
F=inv(C*inv(-A+B*K)*B)
sys_re_filt=ss(A-B*K,B*F,C,D); %sistem dupa F
figure,step(sys_re_filt)
[num,den]=ss2tf(A-B*K,B*F,C,D);
Ho=zpk(tf(num,den))
Hd=Ho/(1-Ho)
Hr=feedback(1,Hd)

Hr=1-Ho
figure,step(100*F*Hr)
%%
syms a b c
Ks=[0 0 1]*inv(ctrb(A,B))*(A-a*eye(3))*(A-b*eye(3))*(A-c*eye(3))
Fs=inv(C*inv(-A+B*Ks)*B)
%%
L=acker(A', C', [-899.132, -899.132, -899.132])'
%Aco  = [-Ra/La, -Ke/La, 0, -K(1)*KPWM/La,-K(2)*KPWM/La,-K(3)*KPWM/La;
        %Kt/J , -Bf/J , 0, 0, 0, 0;
        %0    , 1     , 0, 0, 0, 0;
        %0    , 0     , L(1) -Ra/La, -Ke/La, -L(1);
        %0    , 0     , L(2) Kt/J , -Bf/J , -L(2);
        %0    , 0     , L(3) 0    , 1     , -L(3)]

Aco=[A, -B*K; L*C A-L*C]
Bco  = [F*B; zeros(size(B))];
Cco  = [C, zeros(size(C))];
Dco  = 0;
sys = ss(Aco, Bco, Cco, Dco)
step(sys)
%%
Te=4/3*1e-4;
Hd1=tf([KPWM/La KPWM/La*Bf/J],[1 Ra/La+Bf/J Ra/La*Bf/J+Ke/La*Kt/J])
Hd2=tf([KPWM/La*Kt/J],[1 Ra/La+Bf/J Ra/La*Bf/J+Ke/La*Kt/J])
Hd3=tf([Kt/J*KPWM/La],[1 Ra/La+Bf/J Ra/La*Bf/J+Ke/La*Kt/J 0])
Hr=feedback(1,Hd1*K(1)+Hd2*K(2)+Hd3*K(3))
Hr=1/(1+Hd1*K(1)+Hd2*K(2)+Hd3*K(3))
step(100*F*Hr)
title("Poli -22")

