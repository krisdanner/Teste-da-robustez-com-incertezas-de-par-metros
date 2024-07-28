%|==========================================================|
%| Teste da robustez com incertezas de parâmetros do modelo |
%|==========================================================|
%| Autor: Christian Danner Ramos de Carvalho                |
%|==========================================================|

close all;
clear;
clc;

% Define variáveis e parâmetros
syms h v u m1 Cd rho1

% Outros parâmteros
% m1 = 68.3; % massa inicial kg
m2 = 38.1; % massa final kg
% Cd = 0.43; % coeficiente de arrasto
d = 0.13; % diâmetro do foguete em metros
S = pi*(d^2)/4; % área da seção frontal m^2
% rho1 = 1.1; % densidade do ar (superfície)
rho2 = 0.6; % densidade do ar (3700 metros acima)
T = 25748; % empuxo máximo
g = 9.81; % m/s^2
c1 = 0.5*Cd*rho1*S; % constante aerodinâmica (superfície)
c2 = 0.5*Cd*rho2*S; % constante aerodinâmica (3700 metros acima)

theta = 8; % ângulo de elevação (pitch angle)
phi = 0; % ângulo de rotação (roll angle)

% Define as equações
f1 = v;
f2 = (T*cosd(theta)*u-c1*v^2*sign(v)-m1*g)/m1;

% Linearize a equação em torno do ponto de equilíbrio (0,1080,1)

h_eq = 0;
v_eq = 1080; % velocidade em t = 5s
u_eq = 1;

A = jacobian([f1; f2], [h; v]);
A = subs(A, [h, v, u,], [h_eq, v_eq, u_eq]);
disp(A);

B = jacobian([f1; f2], u);
B = subs(B, [h, v, u], [h_eq, v_eq, u_eq]);
disp(B);

% Parâmetros com incertezas
m1 = ureal('m1',1,'percent',30); % kg (massa inicial com variação de 30%)
Cd = ureal('Cd',1,'percent',10); % coeficiente de arrasto com variação de 10%
rho1 = ureal('rho1',1,'percent',20); % densidade do ar com variação de 20%

A = [0 1;
    0 -(4563*pi*Cd*rho1)/(1000*m1)];

B = [0;
    7008678056137771/(274877906944*m1)];

C = [1 0;
    0 1];

D = 0;

sys = ss(A,B,C,D);
disp(sys);

% Controlador PID
Kp = 24.2; Ki = 280; Kd = 0.523;
s = tf('s');
K = Kp+Ki*1/s+Kd*s;

% Sistema incerto em malha fechada
sysMF = feedback(sys*K, [1 0]); 

% Checar a robustez da margem de estabilidade para o sistema
opt = robOptions('Display','on','Sensitivity','on');
[StabilityMargin, wcu] = robstab(sysMF,opt);

% Menor combinação de parâmetros que desestabiliza o sistema
display(wcu);

wcu.Cd = wcu.Cd*0.95;
sysw = usubs(sysMF,wcu);

% Resposta ao degrau para o sistema com incerteza
figure;
step(sysMF)
title('Resposta ao degrau em MF');
grid;
set(findall(gcf,'type','line'),'linewidth',1);


