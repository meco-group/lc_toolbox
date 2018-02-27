function dx = myODE(t,x,sys,p,u)
%MYODE Summary of this function goes here
%   Detailed explanation goes here

pt = p(t);
ut = u(t);
syspt = sys.f({pt});
syspt = std(syspt{1});

dx = syspt.A*x + syspt.B*ut;
end

