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
figure, step(feedback(Hvth,1),feedback(Hvths,1)); 
legend('Sistem initial','Neglijarea constantei electrice','Location','southeast')
%% întârziere de faza
hold on
nr=50;
z=logspace(-7,1,nr);

for i=1:nr
    for j=1:nr
        Hc=tf([z(i) 1], [z(j) 1]);
        Ho=feedback(series(Hc, Hvths),1);
        [num, den]=tfdata(Ho, 'v');
        x=roots(den);
        y=roots(num);
        ok=1;
        for k=1:length(x)
            if(~isreal(x(k)))
                ok=0;
            end
        end
        if ok
             if(min(abs(x))<(abs(y)))
                plot(z(i), z(j), 'xg')
             elseif(min(abs(x))<=1.000001*(abs(y)))
                plot(z(i), z(j), 'xm')
             elseif(min(abs(x))<=1.000005*(abs(y)))
                plot(z(i), z(j), 'xk')
             elseif(min(abs(x))<=1.00001*(abs(y)))
                plot(z(i), z(j), 'xr')
             elseif(min(abs(x))<=1.0005*(abs(y)))
                plot(z(i), z(j), 'xy')
             elseif(min(abs(x))<=1.0001*(abs(y)))
                plot(z(i), z(j), 'xb')
             end
       end
    end
end
xlabel('Tz') 
ylabel('Tp') 
hold on

x=0.011514;
y=0.000828643;
a=x*y
zz=logspace(-2,1,1000);
zy=zeros(1,1000);

for i=1:1000
    zy(i)=a/zz(i);
end
plot(zz,zy,'m')
%% timp de raspuns
hold on;
for i=1:1000
    Hc=tf([zz(i) 1], [a/zz(i) 1]);
    Ho=feedback(series(Hc, Hvths),1);
    [num, den]=tfdata(Ho, 'v');
    plot(i,min(abs(roots(den))),'xb');
end
xlabel('Tz') 
ylabel('Pol dominant') 

%% simulare
Tz=30*sqrt(a/30)
Hc1=tf([Tz 1], [a/Tz 1])
Ho1=feedback(series(Hvths,Hc1),1)
zpk(Ho1)
figure, step(Ho1,1)
title('Regulator cu avans de fază');
%% nichols
nichols(Hvths,series(Hvths,Hc1))
legend('Fara regulator', 'Cu regulator')
%% senzitivitate
L=series(Hvths,Hc1)
S=1/(1+L);
T=L/(1+L);
KS=Hc1/(1+L);
figure, 
subplot(1,3,1),bode(S)
title('Senzitivitate');
subplot(1,3,2), bode(T)
title('Senzitivitate complementară');
subplot(1,3,3), bode(KS)
title('Efort de control');
%% maxim 5V motor
Tz=30*sqrt(a/30)
Hc1=tf([Tz 1], [a/Tz 1])
Ho1=feedback(Hc1,Hvths);
figure, step(100*Ho1)
title('Intrare motor [V]');
Ho1=feedback(series(Hvths,Hc1),1)
zpk(Ho1)
figure, step(Ho1,10)
title('Regulator cu întârziere de fază');
%% întârziere de faza
hold on
nr=100;
z=logspace(1,3,nr);

for i=1:nr
        Hc=tf([z(i) 1], [z(i)/0.045 1]);
        Ho=feedback(series(Hc, Hvths),1);
        [num, den]=tfdata(Ho, 'v');
        x=roots(den);
        y=roots(num);
        ok=1;
        for k=1:length(x)
            if(~isreal(x(k)))
                ok=0;
            end
        end
        if ok
            plot(z(i),min(abs(x))/(abs(y)),'xb');   
        end
end
xlabel('Tz') 
ylabel('Raport pol dominant-zero') 
%% sistem original
Tz=30*sqrt(a/30)
Hc1=tf([Tz 1], [a/Tz 1])
Ho1=feedback(Hc1,Hvth);
figure, step(100*Ho1)
title('Intrare motor [V]');
Ho1=feedback(series(Hvth,Hc1),1)
zpk(Ho1)
figure, step(Ho1)
title('Regulator cu întârziere de fază');
%% C
Ra_ = ureal('Ra_',Ra,'Percentage',10);
La_ = ureal('La_',La,'Percentage',10);
Ke_ = ureal('Ke_',Ke,'Percentage',10);
Kt_ = ureal('Kt_',Kt,'Percentage',10);
J_ = ureal('J_',J,'Percentage',10);
Bf_ = ureal('Bf_',Bf,'Percentage',10);
KPWM_ = ureal('KPWM_',KPWM,'Percentage',10);
%pentru sistemul simplificat
%numitor H= det(A)
A=[-Ra/La, -Ke/La, 0; 
    Kt/J , -Bf/J , 0; 
    0    , 1     , 0];
B=[KPWM/La; 0; 0 ];
C=[0 0 1];
D=[0 0; 0 0];
%numaratorul l-am calculat cu matricele B si C pentru 2 intrari, respectiv iesiri
%numarator H = C*(s*I3-A')B
H=tf(KPWM_/La_*Kt_/J_-1/J_ , [1 Ra_/La_+Bf_/J_ Ra_/La_*Bf_/J_+Ke_/La_*Kt_/J_ 0]);
Ho1=feedback(series(H,Hc1),1)
figure, step(Ho1,Ho1.NominalValue);
title("Sistem A")