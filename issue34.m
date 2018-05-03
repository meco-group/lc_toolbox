%% reproduce issue 34

load issue34;

G = ss(A,B,C,D);
G = IOSystem(fromstd(G));
K = IOSystem(1,1);
r = Signal(1);
y = G.out();
u = G.in();
e = r - y;
conn = [K.in() == e; K.out() == u];
CL = IOSystem(G,K,conn);

S = Channel(e/r, 'S');
T = Channel(y/r, 'T');

MS = Weight.DC(6);
WS = Weight.LF(1,1,-60);
WT = Weight.HF(3,1,-60);

obj = WS*S;
cstr = [MS*S <= 1; WT*T <= 1];
opts.controller_name = 'K1';
[~,~,info1] = CL.solve(obj,cstr,K,opts);

obj = WS*S;
cstr = [[MS*S ; WT*T] <= 1];
opts.controller_name = 'K2';
[~,~,info2] = CL.solve(obj,cstr,K,opts);

showall(info1,info2);