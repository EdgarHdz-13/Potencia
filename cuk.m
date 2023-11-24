clear 
clc
close all
opts = bodeoptions;
opts.FreqUnits = 'Hz';
s = tf('s');

% Analisis Controlador%
syms VL1 IL1 VC1 IC1 VL2 IL2 VC2 IC2 D RON R L1 L2 C1 C2
Dp = 1-D;
IO = 3;
VO = 5;
VC2 = VO;
vg = 15;
RL1 = 0.260;
RC1 = 1.5;
RL2 = RL1;
RC2 = RC1;
vd = 0.7;
RON = 0.003;
IC2 = 0;
T = 1/50000;

eq13 = R*IO == VO;
R = solve(eq13, R);

eq15 = IL2-VC2/R == 0;
IL2 = solve(eq15,IL2);
IC1 =-IL2;

eq16 = -IL2*D+IL1*Dp == 0;
IL1 = solve(eq16,IL1);

eq17 = vg-RL1*IL1-(RON*(IL1-IC1))*D+(-VC1-RC1*IC1-vd)*Dp == 0;
eq18 = -VC2-RC2*IC2-RL2*IL2+(-RON*(IL1-IC1)+VC1+RC1*IC1)*D-vd*Dp == 0;

[D,VC1] = solve([eq17,eq18],[D,VC1]);
D = D(2);
VC1 = VC1(2);
% Conversión los valores están como 1x1 sym 
syms IL1

R = double(R);
IL2 = double(IL2);
IC1 = double(IC1);
D = double(D);
Dp = (D-1);
VC1 = double(VC1);

eq16 = -IL2*D+IL1*Dp == 0;
IL1 = double(solve(eq16,IL1));


eq1 = VL1 == vg-RL1*IL1-RON*(IL1-IC1);
eq2 = VL2 == -RON*(IL1-IC1)+VC1+RC1*IC1-RL2*IL2-VC2-RC2*IC2;
eq3 = IC1 == IL2;
eq4 = IC2 == IL2-VC2/R;

VL1 = double(solve(eq1,VL1));
VL2 = double(solve(eq2,VL2));
% IC1 == IL2;               % Se cumple
% IC2 = solve(eq4,IC2);     % También se cumple
%
DeIL1 = 0.30;
DeIL2 = 0.30;
DeVC1 = 0.03;
DeVC2 = 0.03;


eq22 = (2*DeIL1)/(D*T) == VL1/L1;
L1 = solve(eq22,L1);

eq23 = (2*DeIL2)/(D*T) == VL2/L2;
L2 = solve(eq23,L2);

eq24 = 2*DeVC1 == (DeIL1*T)/(8*C1);
C1 = solve(eq24,C1);

eq25 = 2*DeVC2/(D*T) == VC2/(R*C2);
C2 = solve(eq25,C2);

% eq7 = (2*0.30)/(0.302*(1/50000)) == (-0.003*(0.695-3)+15.460+1.5*-3-0.260*3-5)/L2;
% L2 = solve(eq7,L2);
% 
% eq8 = 2*0.03 == (0.30*T)/(8*C1);
% C1 = solve(eq8,C1);
% 
% eq9 = (2*0.03)/(D*T)==5/(R*C2);
% C2 = solve(eq9,C2);

% Conversión de valores a double
L1 = double(L1);
L2 = double(L2);
C1 = double(C1);
C2 = double(C2);

% Matrices
Am = [  -((RL1+RON*D)/L1)   ,-(RON*D)/L1        ,(1-D)/L1   ,0;
        -(RON*D)/L2         ,(RL2+RON*D)/L2     ,D/L2       ,-1/L2;
        ((1-D)*T)/C1        ,-D/C1              ,0          ,0;
        0                   ,1/C2               ,0          -1/(C2*R)
    ];

Bm = [  1/L1                ,-((IL1-IL2)*RON)/L1;
        0                   ,-(((IL1-IL2)*RON)+D)/L2;
        0                   ,(-IL2-IL1*T)/C1;
        0                   ,0
     ];
Cm = [0 0 0 1];
Dm = [0 0];

sys = ss(Am,Bm,Cm,Dm);

% Funciones de transferencia
H = tf(sys);
Hvg = H(1); 
HD = H(2); 
% sisotool(HD);
