function [dx] = pendulum_model(t,x,u,m,l,g,J,c)
%PENDULUM_MODEL Summary of this function goes here
%   x = [theta,dtheta,x,dx]

u = u(t);

dx(1,1) = x(2);
dx(2,1) = (m*g*l*sin(x(1)) - c*x(2) + m*l*cos(x(1))*u - m*l*sin(x(1))*x(4)*x(2))/(J+m*l^2);
dx(3,1) = x(4);
dx(4,1) = u;

end

