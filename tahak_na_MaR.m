clc;clear all;close all

%%
s = tf('s')

%% Přenosy (porucha, řízení, regulační odchylka)

G = 5/(s+2)                   % regulovaná soustava
Gr = pid(2,3,0)               % regulátor
Go = series(G,Gr)

Ge = feedback(1,Go)           % přenos regulační odchylky
Gv = feedback(G,Gr)           % přenos poruchy
Gw = feedback(Go,1)           % přenos řízení

% Zápis:

% GE(s) = E(s)/W(s) = 1/(1+Go(s)) = 1/(1+G(s)*GR(s))
% GV(s) = Y(s)/V(s) = G(s)/(1+Go(s)) = G(s)/(1+G(s)*GR(s))
% GW(s) = Y(s)/W(s) = Go(s)/(1+Go(s)) = (G(s)*GR(s))/(1+G(s)*GR(s))

%% Odezva systému se setrvačností (obrácená exponenciála, zesílení)

% y(inf) = hodnota na ose Y, kde se signál ustálí
% u(inf) = zadáno (např. u(t) = 3n(t))

% T = hodnota osy X v 63,2% osy Y
% K = y/u

% Gap = K/(T*s+1)
% G = (K/T)/(T*s+1)


%% Odezva systému (kmitavý průběh,peaky)

% y(inf) = hodnota na ose Y, kde se signál ustálí
% u(inf) = zadáno

% TA = perioda mezi A1 a A2
% A1 = výška 1. peaku od ustálené hodnoty y(inf)
% A2 = výška 2. peaku od ustálené hodnoty y(inf)

% Zbytek dosazuju do vzorců


%% Diferenciální rovnice zapsat jako přenos (laplace,časové konstanty)

% Znám rovnici složenou z y(t) a u(t) - jako zadání Laplaceovy transformace:

% Převedu pomocí záměny y´´(t) => s^2Y(s), y´(t) => sY(s), y(t) => Y(s),
% u´(t) => sU(s), u(t) => U(s)
% Potom si vytknu Y(s) a U(s) -> fám do zlomku G(s) = Y(s)/U(s) a doplním
% hodnoty ze závorek zbylých po vytýkání -> hodnoty potom dodám do matlabu

cit = [6 1]                   % čitatel (u(t))
jmen = [30 11 1]              % jmenovatel (y(t))
G = tf(cit,jmen)
[Z,P,K] = tf2zp(cit,jmen)

Gzpk = zpk(Z,P,K)             % tady výsledek pro zápis pomocí rozložení nul a pólů

% Zapíšu si Gzpk a z něj určuji časové konstanty

% určování konstant:
                              % tady výsledek pro zápis pomocí časových konstant
% vyjde přenos v určitém tvaru -> v závorkách musí být na pravé straně
% závorky vždy "+1" (musí tam být celá 1) -> upravím tak, že to číslo na
% místě, kde by měla být celá 1, vytknu před závorku -> nejspíš se mi to
% pak vykrátí -> potom na čísla, která jsou násobena konstantou "s",
% aplikuju vzorce pro T = |1/p| a Tau = |1/n| -> dostanu výsledné hodnoty
% nul a pólů

%% Póly a nuly

% p = -(1/T)
% n = -(1/Tau)


%% Převod PID na PSD (časová konstanta)

G = (5*s+1)/((s+10)*(s+5))    % regulovaná soustava
Gr = pid(2,0.4,5)             % regulátor

% Vyberu časovou konstantu Ts:

% Určím tak, že vytáhnu póly z regulované soustavy 
% pomocí vzorce p = 1/s a projedu s ním všechny póly regulované soustavy
% -> dostanu výsledky v záporných hodnotách -> potom výsledné záporné
% hodnoty projedu vzorcem T = |1/p| -> nejmenší/největší hodnotu
% vynásobím/podělím 10

% Výsledky:
% p1 = -10, p2 = -5 --> T1 = 0.1, T2 = 0.2 --> vyberu tu nejmenší a udělám
% ji 10x menší -> Ts = 0.01

Ts = 0.01

% Na tohle nesahám
Td = Gr.Kd/Gr.Kp              % derivařní časová konstanta regulátoru
Kp = Gr.Kp                    % proporcionální konstanta regulátoru
Ti = Gr.Kp/Gr.Ki              % integrační časová konstanta regulátoru
d0 = Kp*(1+Ts/Ti+Td/Ts)
d1 = -Kp*(1+2*Td/Ts)
d2 = Kp*Td/Ts

% přírůstkový PSD regulátor
Grz = tf([d0 d1 d2],[1 -1 0],Ts,'variable','z^-1')

% Úprava na diferenční rovnici (znaménka podle výsledku):

% p1*U(z)-p2z^-1*U(z) = n1*E(z) - n2z^-1*E(z) + n3z^-2*E(z)  % n = nuly -> hodnoty z čitatele, p = póly -> hodnoty ze jmenovatele

% u[k] = u[k-1] + n1*e[k] - n2*e[k-1] + n3*e[k-2]            % tohle je výsledek


%% Kvalita regulace v časové oblasti (přenos s překmitem, dobou náběhu atd)

G = (5*s+1)/((s+10)*(s+5))    % regulovaná soustava
Gr = pid(2,0.4,5)             % regulátor

Go = series(G,Gr)             % přenos soustavy
Gw = feedback(Go,1)           % přenos regulované veličiny (řízení)
figure,step(Gw)               % step response
stepinfo(Gw)                  % Informace jako rise time, settling time, overshoot atd.


%% integrální kritéria + kvalita regulace v kmitočtové oblasti (kvadratická, ITAE)

G = 2/((6*s+1)*(5*s+1))                 % regulovaná soustava
Gr = ((6*s+1)*(5*s+1))/(s*(0.1*s+1))    % regulátor
Go = series(G,Gr);
Ge = feedback(1,Go)

t_stop = 3;                   % dostatečně dlouhý čas -> volím tak, aby bylo vidět, že se signál ustálil

figure(),step(Ge,t_stop)      % pro ověření "t_stop"

[e,t] = step(Ge,t_stop);

JK = trapz(t,e.^2)            % kvadratické kritérium (integrál od 0 do inf z e^2(t) dt)
JITAE = trapz(t,t.*abs(e))    % ITAE kritérium (integrál od 0 do inf z t*|e|(t) dt)

% kvalita v kmitočtové oblasti
figure(),margin(Go)
[Gm,pm,omega_faze,omega_rezu] = margin(Go)  % Gm = amplitudová bezpečnost, Pm = fázová bezpečnost
