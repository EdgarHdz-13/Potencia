clear 
clc
close all
opts = bodeoptions;
opts.FreqUnits = 'Hz';
s = tf('s');

% Analisis Controlador%
syms VL1 IL1 VC1 IC1 VL2 IL2 VC2 IC2 R D RON RO IO VO L1 L2 C1 C2
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

eq1 = RO*IO == VO;
RO = solve(eq1, RO);

eq2 = IL2-VC2/RO == 0;
IL2 = solve(eq2,IL2);
IC1 =-IL2;

eq3 = -3*D+IL1*Dp == 0;
IL1 = solve(eq3,IL1);

eq4 = vg-RL1*IL1-(RON*(IL1-IC1))*D+(-VC1-RC1*IC1-vd)*Dp == 0;
eq5 = VC2-RC2*IC2-RL2*IL2+(-RON*(IL1-IC1)+VC1+RC1*IC1)*D-vd*Dp == 0;

[D,VC1] = solve([eq4,eq5],[D,VC1]);
D = D(2);
VC1 = VC1(2);
D = 0.302;
eq6 = (2*0.30)/(0.302*(1/50000)) == (15-RL1*0.695-0.003*(0.695-3))/L1;
L1 = solve(eq6,L1);

eq7 = (2*0.30)/(0.302*(1/50000)) == (-0.003*(0.695-3)+15.460+1.5*-3-0.260*3-5)/L2;
L2 = solve(eq7,L2);

eq8 = 2*0.03 == (0.30*T)/(8*C1);
C1 = solve(eq8,C1)

eq9 = (2*0.03)/(D*T)==5/(RO*C2);
C2 = solve(eq9,C2)

