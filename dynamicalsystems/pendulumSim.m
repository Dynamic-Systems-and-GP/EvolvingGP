function [y,x1,u] = pendulumSim(x0,u,Ts,params)
    Ps=params;
    f=@(t,x) odefun(t,x,u,Ps);
    x0ode=[x0.z x0.dz x0.fi x0.dfi]';
    [t,x1ode]=ode23(f,[0 Ts],x0ode); %11 steps per each run is minimum
    x1.z =x1ode(end,1);
    x1.dz=x1ode(end,2);
    x1.fi =x1ode(end,3);
        
    %output:
    y.z=x1.z;
    y.fi=x1.fi;
end



function [dx]=odefun(t,x,u,Ps)
dx=NaN(3,1);
m=Ps.m;
l=Ps.l;
M=Ps.M;
g=Ps.g;

z =x(1);
dz=x(2);
fi =x(3);
dfi =x(4);

dx(1)=dz;
dx(2)=(- l*m*sin(fi)*dfi^2 + F + g*m*cos(fi)*sin(fi))/(M + m - m*cos(fi)^2);
dx(3)=dfi;
dx(4)=(- l*m*cos(fi)*sin(fi)*dfi^2 + F*cos(fi) + g*m*sin(fi) + M*g*sin(fi))/(l*(M + m - m*cos(fi)^2));

end

%   rez=solve(...
%   '(M+m)*ddx-m*l*ddfi*cos(fi)+m*l*dfi^2*sin(fi)==F',...
%   'l*ddfi-g*sin(fi)=ddx*cos(fi)',...
%   ddx,ddfi)