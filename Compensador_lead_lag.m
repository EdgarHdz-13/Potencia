%% Principal
clear 
close all
opts = bodeoptions;
opts.FreqUnits = 'Hz';
s = tf('s');

%% Solve
syms R iL2 D vc1 iL1 L1 L2 C1 C2 Ron RL1 RL2 T vc2 vg ig vd
syms ic1
Dp = 1-D;

eq13    = R*ig == vg;
eq12    = iL2-vc2/R == 0;
eq11    = -iL2*D+iL1*(1-D) == 0;
eq9     = vg - RL1*iL1-((iL1-iL2)*Ron)*D+(-vc1-vd)*Dp == 0;
eq10    = -vc2-RL2*iL2+(vc1-(iL1-iL2)*Ron)*D-vd*Dp == 0;    

% Ecuaciones
R = solve(eq13, R);                
iL2 = solve(eq12, iL2);       
iL1 = solve(eq11, iL1);
[D,vc1] = solve([eq9,eq10],[D,vc1]);

%% Resolver las ecuaciones
% Variables
vg = 5; ig = 3; vc2 = 5; Ron = 30e-3;

% Ecuaciones de nuevo
syms R iL2 iL1 D vc1

eq13    = R*ig == vg;
R = solve(eq13, R);  

eq12    = iL2-vc2/R == 0;
iL2 = solve(eq12, iL2); 

eq11    = -iL2*D+iL1*(1-D) == 0;
iL1 = solve(eq11, iL1);

eq9     = vg - RL1*iL1-((iL1-iL2)*Ron)*D+(-vc1-vd)*Dp == 0;
eq10    = -vc2-RL2*iL2+(vc1-(iL1-iL2)*Ron)*D-vd*Dp == 0;    
[D,vc1] = solve([eq9,eq10],[D,vc1]); 


%% Valores del sistema
R   =   5/3; 
iL2 =   5/R; 
D   =   0.2786;
vc1 =   20.0558;
iL1 =   (iL2*D)/(1-D);
L1  =   209e-6;
L2  =   209e-6;
C1  =   25e-6;
C2  =   836e-6;
Ron =   13.9e-3;
RL1 =   30e-3;
RL2 =   30e-3;
T   =   1/50e3;

%% Matrices
Am = [  -((RL1+Ron*D)/L1)   ,-(Ron*D)/L1        ,(1-D)/L1   ,0;
        -(Ron*D)/L2         ,(RL2+Ron*D)/L2     ,D/L2       ,-1/L2;
        ((1-D)*T)/C1        ,-D/C1              ,0          ,0;
        0                   ,1/C2               ,0          -1/(C2*R)
    ];

Bm = [  1/L1                ,-((iL1-iL2)*Ron)/L1;
        0                   ,-(((iL1-iL2)*Ron)+D)/L2;
        0                   ,(-iL2-iL1*T)/C1;
        0                   ,0
     ];
Cm = [0 0 0 1];
Dm = [0 0];

sys = ss(Am,Bm,Cm,Dm);

%% Funciones de transferencia
H = tf(sys);
H1 = H(1); 
H2 = H(2); 
H3 = H(3); 

% sisotool(H1);

% bode(H1, opts); title('H1(s)');
% figure;
% bode(H2, opts); title('H2(s)');
% figure;
% bode(H3, opts); title('H3(s)');

%% Diseño del compensador lead

fdis    = 4e3;         %Frecuencia deseada
fo      = 631;         %Frecuencia de resonancia
phase    = 60;    
theta   = (phase*pi)/180;

fp = fdis*(sqrt((1+sin(theta))/(1-sin(theta))));
fz = fdis*(sqrt((1-sin(theta))/(1+sin(theta))));

gain = 25.53;          %Ganancia del sistema
K = ((fdis/fo)^2)*(1/gain)*sqrt(fz/fp);

%% Función de transferencia del LEAD
wz = 2*pi*fz;
wp = 2*pi*fp;

CsLead = K*((s/wz)+1)/((s/wp)+1);

% Función de transferencia deseada
TsLead = -CsLead*H1;

% Plot de la función de transferencia
figure;
bode(TsLead, opts);hold on;
bode(H1, opts, 'K--');  hold off;

%% Diseño del compensador lag
fc = 6e3;
fL = 10;

wLLag  = 2*pi*fL;
wcLag = 2*pi*fc;
woLag = 2*pi*fdis;

Gco = (((wcLag/woLag)^2)/K);
CsLag = Gco*(1 + (wLLag/s));

% Función de transferencia deseada
TsLag = CsLag*H1;
figure;
bode(Hvd, opts, 'K--'); hold on;
bode(TsLag, opts); hold off;

%% Lead-lag

CsLeadLag = K*( ( (s/wz)+1)*((wLLag/s) + 1 )/ ((s/wp)+1) );

TsLeadLag = CsLeadLag*Hvd;
bode(Hvd, opts,'K--'); hold on;
bode(TsLeadLag, opts); hold off;

