function y = myOutput(t,x,sys,p,u)
%MYOUTPUT Summary of this function goes here
%   Detailed explanation goes here

pt = p(t);
ut = u(t);

syspt = sys.f({pt});
syspt = std(syspt{1});

y = syspt.C*x + syspt.D*ut;
end

