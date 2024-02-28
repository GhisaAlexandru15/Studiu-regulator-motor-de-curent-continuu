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
B=[KPWM/La 0; 0 -1/J; 0 0];
C=[0 1 0; 0 0 1];
D=[0 0; 0 0];
H=tf(ss(A,B,C,D))
figure,step(H)
Hvth=zpk(minreal(H(2,1)))
[num, den]=tfdata(Hvth, 'v');
rden=roots(den);
Hvths=zpk(tf(num(4)/(-rden(2)), [1 -rden(3) 0]))
%
Tz=500
Hc1=tf([Tz 1], [Tz/0.045 1]);
Te=4/3*1e-4;
Hcd=c2d(Hc1,Te,'zoh')
[num,den]=tfdata(Hcd,'values')
Hcd10=tf([0.0450000000  -0.0449999880],[1.0000000000  -0.9999999880],Te)
Hcd7=tf([0.0450000  -0.0449999],[1.0000000  -0.9999999],Te)
Hcd5=tf([0.04500  -0.04499],[1.00000  -0.99999],Te)
Hcd3=tf([0.045 -0.044],[1.000  -0.999],Te)

Hfd=c2d(Hvths,Te,'zoh')
Hod=feedback(series(Hcd, Hfd),1);
Hod10=feedback(series(Hcd10, Hfd),1);
Hod7=feedback(series(Hcd7, Hfd),1);
Hod5=feedback(series(Hcd5, Hfd),1);
Hod3=feedback(series(Hcd3, Hfd),1);

step(Hod,Hod10,Hod7,Hod5,Hod3)
legend("15 zecimale","10 zecimale","7 zecimale","5 zecimale","3 zecimale")
%%
Hvthd=c2d(Hvth,Te,'zoh');
Ho1=feedback(Hcd7,Hvthd);
figure, step(100*Ho1)
title('Intrare motor [V]');
Ho1=feedback(series(Hvthd,Hcd7),1)
zpk(Ho1)
figure, step(Ho1,Hod7)
title('Regulator cu întârziere de fază');