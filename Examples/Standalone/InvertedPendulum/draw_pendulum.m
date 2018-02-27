function [ax] = draw_pendulum(ax,x,theta,varargin)
%PENDULUM_MOVIE Summary of this function goes here
%   Detailed explanation goes here    

if(nargin>3)
    l = varargin{1}.l;
    wc = varargin{1}.wc;
    hc = varargin{1}.hc;
    wp = varargin{1}.wp;
    hp = varargin{1}.hp;
end

y0 = 0;
r = hc/4;
N = 100;

c = linspace(0,2*pi,N)';
Pwheel1 = repmat([x-l/2 y0+r],[N,1]) + [cos(c) sin(c)]*r;
Pwheel2 = repmat([x+l/2 y0+r],[N,1]) + [cos(c) sin(c)]*r;
Pcart = [x-wc/2 y0+r;...
         x-wc/2 y0+hc+r;...
         x+wc/2 y0+hc+r;...
         x+wc/2 y0+r;...
         x-wc/2 y0+r];
     
a = linspace(0,pi,N)';     
Ppendle_upper = repmat([0 hp-wp+r+hc/2],[N,1]) + [-cos(a) sin(a)]*wp/2;
Ppendle_lower = repmat([0 r+hc/2],[N,1]) + [cos(a) -sin(a)]*wp/2;
Ppendle = [-wp/2,r+hc/2;...
           -wp/2,hp-wp+r+hc/2;...
           Ppendle_upper;...
           wp/2,hp-wp+r+hc/2;...
           wp/2,r+hc/2;...
           Ppendle_lower;...
           -wp/2,r+hc/2];
       
Rot=[cos(theta) sin(theta);-sin(theta) cos(theta)];
Ppendle = Ppendle*Rot + repmat([x y0],[5+2*N,1]);

hold(ax,'on');
% plot(ax,Pcart(:,1),Pcart(:,2),'k-','Linewidth',2);
% plot(ax,Pwheel1(:,1),Pwheel1(:,2),'k-','Linewidth',2);
% plot(ax,Pwheel2(:,1),Pwheel2(:,2),'k-','Linewidth',2);
%plot(ax,Ppendle(:,1),Ppendle(:,2),'b-','Linewidth',2);
patch(Pcart(:,1),Pcart(:,2),[139 137 137]/255,'EdgeColor','k','Linewidth',2)
patch(Pwheel1(:,1),Pwheel1(:,2),[139 137 137]/255,'EdgeColor','k','Linewidth',2)
patch(Pwheel2(:,1),Pwheel2(:,2),[139 137 137]/255,'EdgeColor','k','Linewidth',2)
patch(Ppendle(:,1),Ppendle(:,2),[65 105 255]/255,'EdgeColor','b','Linewidth',2);
hold(ax,'off');

axis(ax,'equal');
if(nargin>4)
    axis(ax,varargin{2});
end
% daspect(ax,[1 1 1]);

end

