%% Funzioni di inizializzazione
clear      % cancella le variabili nel workspace
close all  % chiude tutte le figure
clc        % pulisce la riga di comando
%% Impostazioni per le figure
font_size   = 15;   % dimensione dei font
marker_size = 12;   % dimensione dei marker
line_width  = 2;    % dimensione delle linee
plot_x0     = 500;  % coordinata x delle figure
plot_y0     = 300;  % coordinata y delle figure
plot_width  = 600;  % larghezza delle figure
plot_height = 400;  % altezza delle figure
%% Variabile per definire funzioni di trasferimento
s = tf('s');
%% Impostazioni per Simulink
% Parametri del solver: Modeling --> Model Settings --> Solver
MaxStep  = 1e-5; % massimo step di simulazione
                 % in genere, prendere questo parametro almeno 10 volte più
                 % piccolo della costante di tempo più veloce da simulare
RelTol   = 1e-5; % tolleranza relativa del solver
% Tempo di simulazione: Stop Time in alto a destra della finestra
t_f      = 10;   % scegliere un valore tale da catturare tutta la risposta
%% Parametri elettrici
Re = 4.8;
Le = 55e-3;
%% Parametri termici
T0 = 273.15;
c = 6.2;
l = 15e-2;
d = 0.2e-3;
Cth = 25.1e-3;
Rth = 16.7;
Tp = 15+T0;
h = 6.55;
Tamb = 24+T0;
%% Parametri meccanici
r1 = 0.45e-2;
r2 = 4.5e-2;
J = 2e-4;
m = 15e-3;
beta = 1.5;
%% Parametri K(T)
Kmax = 3.92e3;
Tavg = 70;
Tdiff = 20;
x3_eq = 4e-2;
%% Parametri costruiti
m_eq = J/r2 + m*r2;
g = 9.81;
Ramb = 1/(h*pi*d*l);
S = pi*d*l;
rho = r2/r1;
deltal_max = 4/100*l;

%% Calcolo punti di equilibrio
K_star = (m*g*r2)/(r1*(deltal_max-x3_eq*r1/r2));
x2_eq = (Tavg+(Tdiff/c)*log(K_star/(Kmax-K_star)));
x2_eq_k = x2_eq +T0;
x1_eq = sqrt((1/Re)*(h*S*(x2_eq_k-Tamb)+(1/Rth)*(x2_eq_k-Tp)));
x4_eq = 0;
expo = exp(c*((x2_eq_k-(Tavg+T0)))/(Tdiff+T0));
K_diff = (c*Kmax/(Tdiff+T0))*(expo/((1+expo)^2));
K_x2eq = Kmax*(1-(1/(1+expo)));
u_eq = Re*x1_eq;

    % Stampa delle soluzioni
    disp('Soluzione del vettore colonna:');
    disp('Ie* = ');
    disp(x1_eq);
    disp('T* = ');
    disp(x2_eq);
    disp('Z* = ');
    disp(x3_eq);
    disp('z* = ');
    disp(x4_eq);
    disp('V* = ');
    disp(u_eq);

%% Definizione modello linearizzato 
A = [-Re/Le 0 0 0;
    2*Re*x1_eq/Cth -1/Cth*(h*S+1/Rth) 0 0;
    0 0 0 1;
    0 ((r1*(deltal_max-x3_eq*r1/r2))/(r2*m_eq))*K_diff -((r1^2)/((r2^2)*m_eq))*K_x2eq -beta/m_eq];
B = [1/Le; 0; 0; 0];
C = [0 0 1 0];
D = 0;

LTI_model = ss(A, B, C, D);
G = tf(LTI_model);

%% Specifiche statiche e dinamiche
% errore nullo, andremo quindi a inserire un polo nell'origine di L(s);
Ta = 1.7;
% Sovraelongazione percentuale < 5%
%% Costruzione FdT per inner loop (Gt*Ge linearizzate in Ie*)
dAt = [-Re/Le               0;
        2*Re*x1_eq/Cth    -1/(Rth*Cth)-1/(Ramb*Cth)];
dBt = [1/Le   0;
       0     1/(Rth*Cth)];
dCt = [0 1];
dDt = 0;
model_t = ss(dAt, dBt, dCt, dDt);
Gt = tf (model_t);
[nt1, dt1] = tfdata (Gt(1), 'v'); % Gt(1) lega V e T, mentre Gt(2) lega dTp e uscita
Gt(1)

%% Costruzione FdT per outer loop (Gm linearizzata in T* e in x3*)
n = exp (6.2/Tdiff*(x2_eq-Tavg));
dAm = [0                                1;
     -Kmax * (1 - 1 / (1 + n)) * r1/(m_eq*rho)    -beta*r2/m_eq];
dBm = [0;
       Kmax*r1/m_eq*(4/100*l-x3_eq/rho)*6.2/Tdiff*n/(1+n)^2];
dCm = [1 0];
dDm = 0;
model_lin_m = ss(dAm, dBm, dCm, dDm);
Gm = tf (model_lin_m); 
Gm

%% Progetto regolatore inner-loop: PID
wc1 = 290;
Mf1 = 65;
%Costruzione regolatore
[M1,P1,W1] = bode(Gt(1)/(s)); 
[V,i] = min(abs(W1-wc1)); 
ModGt = 20*log10(M1(i));
ArgGt = P1(i);
Pdi = (-180+Mf1-ArgGt); 
rho_star1 = tand(Pdi/2);
tauz = rho_star1/wc1;
A1 = 20*log10(sqrt (1 + rho_star1^2));
u1 = 10^(-(2*A1+ModGt)/20);
taup = 1/(10*wc1);
R1 = u1*(1 + tauz*s)*(1 + tauz*s)/s/(1 + taup*s)/(1+0.000012*s);
L1 = R1*Gt(1);
%Funzioni di sensitività
F1 = minreal(L1/(L1+1));
S1 = minreal (1/(1 + L1));
Q1 = minreal (R1/(1 + L1));
%Modello nello spazio degli stati
mod_R1 = ss(R1);
A1 = mod_R1.A;
B1 = mod_R1.B;
C1 = mod_R1.C;
D1 = mod_R1.D;

%Numeratore e denominatore del regolatore
[n1, d1] = tfdata (R1, 'v');
R1

%% Progetto regolatore outer-loop: PI
Mf = 65;
wc = 460/(Ta*Mf);
%Costruzione del regolatore
[M2,P2,W2] = bode(Gm/s); 
[V,i]=min(abs(W2-wc)); 
ModGm = 20*log10(M2(i));
ArgGm = P2(i);
Pdo = (-180+Mf-ArgGm);
rho_star = tand(Pdo);
tau = rho_star/wc;
Ao = 20*log10(sqrt(1+rho_star^2));
u = 10^(-(Ao+ModGm)/20);
R2 = u*(1+tau*s)/(s*(1+0.024*s)*(1+0.0001*s));
L2 = R2*Gm;
%Funzioni di sensitività
F2 = minreal(R2*F1*Gm/(1 + R2*F1*Gm));
S2 = minreal(1/(1 + R2*F1*Gm));
Q2 = minreal(R2/(1 + R2*F1*Gm));
R2;
%Modello nello spazio degli stati
mod_R2 = ss(R2);
A2 = mod_R2.A;
B2 = mod_R2.B;
C2 = mod_R2.C;
D2 = mod_R2.D;

%% Stampa dei grafici
%Valutazione dello sforzo di controllo dell'inner loop
figure(1)
set(1,'position',[plot_x0,plot_y0,plot_width,plot_height])
bodemag(Gt(1), R1, Q1)
grid on
h = findobj(1, 'Type', 'Line');
set(h, 'LineWidth', line_width)
title('Valutazione della sensitività del controllo (inner loop)', 'FontSize', font_size)
legend('Gm', 'Rm', 'Qm', 'FontSize', font_size)

%Valutazione dello sforzo di controllo dell'outer loop
figure(2)
set(2,'position',[plot_x0,plot_y0,plot_width,plot_height])
bodemag(Gm, R2, Q2)
grid on
h = findobj(2, 'Type', 'Line');
set(h, 'LineWidth', line_width)
title('Valutazione della sensitività del controllo (outer loop)', 'FontSize', font_size)
legend('Gv', 'Rv', 'Qv', 'FontSize', font_size)

%% Diagramma di bode dell'inner loop (L monte)
w = logspace(-1, 4, 1000);

% Calcolo della risposta in frequenza
[mag, phase] = bode(L1, w);
% Creazione della schermata
figure(3);
% Diagramma di Bode del modulo
subplot(2, 1, 1);
semilogx(w, 20*log10(squeeze(mag)));
hold on;

% Disegno del rettangolo rosso 1
w1 = [0.1 0.1 5 5];
mag1 = [60 -60 -60 60];
patch(w1, mag1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

title('Diagramma di Bode (inner loop)');
ylabel('Ampiezza (dB)');
grid on; box on;

% Diagramma di Bode della fase
subplot(2, 1, 2);
semilogx(w, squeeze(phase));
hold on;

xlabel('Frequenza (rad/s)');
ylabel('Fase (gradi)');
grid on; box on;
legend ("L monte")

%% diagramma di bode dell'outer loop
w = logspace(-1, 2, 1000);
% Calcolo della risposta in frequenza
[mag, phase] = bode(L2, w);
% Creazione della schermata
figure(4);

% Diagramma di Bode del modulo
subplot(2, 1, 1);
semilogx(w, 20*log10(squeeze(mag)));
hold on;

% Disegno della zona proibita 1
w1 = [0.1 0.1 4 4];
mag1 = [0 -60 -60 0];
patch(w1, mag1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Disegno della zona proibita 2
w2 = [20 20 100 100];
mag2 = [60 -20 -20 60];
patch(w2, mag2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

title('Diagramma di Bode (outer loop)');
ylabel('Ampiezza (dB)');
grid on;

% Diagramma di Bode della fase
subplot(2, 1, 2);
semilogx(w, squeeze(phase));
hold on;

% Disegno della zona proibita nel diagramma della fase
w3 = [4 4 20 20];
phase3 = [-115 -180 -180 -115];
patch(w3, phase3, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Frequenza (rad/s)');
ylabel('Fase (gradi)');
grid on;
legend ("L valle");

%% Funzioni di sensitività dell'inner loop
figure(5)
set(5,'position',[plot_x0,plot_y0,plot_width,plot_height])
bodemag(F1, S1, R1)
grid on
h = findobj(5, 'Type', 'Line');
set(h, 'LineWidth', line_width)
title('Valutazione delle funzioni di sensitività(inner loop)', 'FontSize', font_size)
legend('F1', 'S1', 'R1', 'FontSize' ,font_size)

%% Funzioni di sensitività dell'outer loop
figure(6)
set(6,'position',[plot_x0,plot_y0,plot_width,plot_height])
bodemag(F2, S2, R2)
grid on
h = findobj(6, 'Type', 'Line');
set(h, 'LineWidth', line_width)
title('Valutazione delle funzioni di sensitività(outer loop)', 'FontSize', font_size)
legend('F2', 'S2', 'R2', 'FontSize' ,font_size)

%% Risposta al gradino (inner loop)
figure(7);
set(7,'position',[plot_x0,plot_y0,plot_width,plot_height]);
[y, t] = step(F1); % salva in y la risposta e in t il tempo
plot(t, y, 'LineWidth', line_width);
grid on;
xlim([0, t(end)]);
xlabel('Tempo [s]', 'FontSize', font_size);
ylabel('Temperatura [K]', 'FontSize', font_size);
title('Risposta al gradino (I-L)', 'FontSize', font_size);

%% Risposta al gradino (outer loop)
figure(8);
set(8,'position',[plot_x0,plot_y0,plot_width,plot_height]);
[y, t] = step(F2); % salva in y la risposta e in t il tempo
plot(t, y, 'LineWidth', line_width);
grid on;
xlim([0, t(end)]);
xlabel('Tempo [s]', 'FontSize', font_size);
ylabel('Posizione [m]', 'FontSize', font_size);
title('Risposta al gradino (O-L)', 'FontSize', font_size);

