clear 
clc
close all
opts = bodeoptions;
opts.FreqUnits = 'Hz';
s = tf('s');

% Variables simbolicas del Jacobiano
syms VL1J IL1J VC1J IC1J VL2J IL2J VC2J IC2J DJ RONJ RJ L1J L2J C1J C2J vgJ RL1J vdJ RL2J
syms TJ
DpJ = 1-DJ;

% Analisis Controlador
syms VL1 IL1 VC1 IC1 VL2 IL2 VC2 IC2 D RON R L1 L2 C1 C2 vg RL1 vd RL2
syms T
Dp = 1-D;

IO = 3;
VO = 5;
VC2 = VO;
vg = 120;
RL1 = 0.260;
RL2 = RL1;
vd = 0.7;
RON = 0.003;
IC2 = 0;
T = 1/50000;

% D. Diseño del convertidor

eq13 = R*IO == VO;                          
R = solve(eq13, R);

eq15 = IL2-VC2/R == 0;                     
IL2 = solve(eq15,IL2);
IC1 =-IL2;

eq16 = -IL2*D+IL1*Dp == 0;                  
IL1 = solve(eq16,IL1);

eq9 = vg-RL1*IL1-(RON*(IL1-IC1))*D+(-VC1-vd)*Dp == 0;
eq10 = -VC2-RL2*IL2+(-RON*(IL1-IC1)+VC1)*D-vd*Dp == 0;

[D,VC1] = solve([eq9,eq10],[D,VC1]);
D = D(2);
VC1 = VC1(2);
% Conversión los valores están como 1x1 sym 
syms IL1

R = double(R);
IL2 = double(IL2);
IC1 = double(IC1);
D = double(D);
Dp = (1-D);
VC1 = double(VC1);

eq16 = -IL2*D+IL1*Dp == 0;                  %REDEFINIR EQ16
IL1 = double(solve(eq16,IL1));


eq1 = VL1 == vg-RL1*IL1-RON*(IL1-IC1);
eq2 = VL2 == -RON*(IL1-IC1)+VC1-RL2*IL2-VC2;
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

% Conversión de valores a double
L1 = double(L1);
L2 = double(L2);
C1 = double(C1);
C2 = double(C2);


% Matrices Normales
% Am = [  -((RL1+RON*D)/L1)   ,-(RON*D)/L1        ,(1-D)/L1   ,0;
%        -(RON*D)/L2         ,(RL2+RON*D)/L2     ,D/L2       ,-1/L2;
%        ((1-D)*T)/C1        ,-D/C1              ,0          ,0;
%        0                   ,1/C2               ,0          -1/(C2*R)
%    ];
% 
% %
% 
% Bm = [  1/L1                ,-((IL1-IL2)*RON)/L1;
%         0                   ,-(((IL1-IL2)*RON)+D)/L2;
%         0                   ,(-IL2-IL1*T)/C1;
%         0                   ,0
%      ];
% Cm = [0 0 0 1];
% Dm = [0 0];

% Jacobiano

% Ecuaciones 
eq28 = (1/L1J)*( vgJ-RL1J*IL1J-(RONJ*(IL1J+IL2J))*DJ+(-VC1J-vdJ)*DpJ );
eq29 = (1/L2J)*( -VC2J-RL2J*IL2J+(-RONJ*(IL1J+IL2J)+VC1J)*DJ-vdJ*DpJ );
eq30 = (1/C1J)*(-IL2J*DJ+IL1J*DpJ);
eq31 = (1/C2J)*( IL2J-VC2J/RJ);

% Función, estados, entrada y salida
f = [eq28,eq29,eq30,eq31];
x = [IL1J, IL2J, VC1J, VC2J];
u = [vgJ, DJ];
y = VC2J;

% Matrices
AmJ = jacobian(f,x);
BmJ = jacobian(f,u);
CmJ = jacobian(y,x);
DmJ = jacobian(y,u);
% Matrices a Sustituir
subsm = [VL1 IL1 VC1 IC1 VL2 IL2 VC2 IC2 D RON R L1 L2 C1 C2 vg RL1 vd  RL2 T];
jsubsm = [VL1J IL1J VC1J IC1J VL2J IL2J VC2J IC2J DJ RONJ RJ L1J L2J C1J C2J vgJ RL1J vdJ RL2J TJ];

Am = double(subs(AmJ, jsubsm,subsm));
Bm = double(subs(BmJ, jsubsm,subsm));
Cm = double(subs(CmJ, jsubsm,subsm));
Dm = double(subs(DmJ, jsubsm,subsm));

sys = ss(Am,Bm,Cm,Dm);

% Funciones de transferencia
H = tf(sys);
Hvc2vg = H(1); 
Hvc2D = H(2); 
% sisotool(Hvc2D);       
% sisotool(Hvc2vg);
% bode(Hvc2D, opts);
% DISEÑO

SysGain = Hvc2D.numerator{1,1}(5)/Hvc2D.denominator{1,1}(5);
freqMax = 30e3;

% %% Diseño del compensador lead
% 
% fcdis           = freqMax/20;           %Frecuencia deseada       1 decada antes
% fo              = 800;                 %Frecuencia de resonancia PSIM
% phasemargin     = 60;                 
% theta           = (phasemargin*pi)/180;
% 
% fp = fcdis*(sqrt((1+sin(theta))/(1-sin(theta))));
% fz = fcdis*(sqrt((1-sin(theta))/(1+sin(theta))));
% 
% K = ((fcdis/fo)^2)*(1/SysGain)*sqrt(fz/fp);
% 
% % Función de transferencia del LEAD
% wz = 2*pi*fz;
% wp = 2*pi*fp;
% 
% CsLead = K*((s/wz)+1)/((s/wp)+1);
% 
% % Función de transferencia deseada
% TsLead = CsLead*Hvc2D;
% 
% % Plot de la función de transferencia
% figure;
% title("Lead Hvc2D");
% bode(TsLead, opts);hold on;
% bode(Hvc2D, opts, 'K--');  hold off;
% 
% %% Diseño del compensador lag
% fc = freqMax/20;
% fL = 15;
% 
% wLLag  = 2*pi*fL;
% wcLag = 2*pi*fc;
% woLag = 2*pi*fcdis;
% 
% Gco = (((wcLag/woLag)^2)/K);
% CsLag = Gco*(1 + (wLLag/s));
% 
% % Función de transferencia deseada
% TsLag = CsLag*Hvc2D;
% figure;
% title("Lag Hvc2D");
% bode(Hvc2D, opts, 'K--'); hold on;
% bode(TsLag, opts); hold off;
% 
% %% Lead-lag
% 
% CsLeadLag = K*( ( (s/wz)+1)*((wLLag/s) + 1 )/ ((s/wp)+1) );
% 
% TsLeadLag = CsLeadLag*Hvc2D;
% 
% figure;
% title("Lead-Lag Hvc2D");
% bode(Hvc2D, opts,'K--'); hold on;
% bode(TsLeadLag, opts); hold off;

% Diseño a mano

% sisotool(Hvc2D); 

CsLeadLagManual = 1.6162*((0.55*s+1)*(0.00012*s+1))/(s*(0.01*s+1));
TsLeadLagManual = CsLeadLagManual*Hvc2D;
bode(TsLeadLagManual, opts); hold off;

