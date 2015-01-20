classdef pendulum_model < process_interface
    %RAKETA_PROCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess=public,SetAccess=public)
        Tsim=110;
        Tper=100;
        %realnocasovne simulacije ne bomo strogo uporabljali s timestamp
        timestamp=NaN;
        %parametri iz modela:
        u=NaN;
        y=NaN;
        x=struct();
        
        b=0.1;  %   damping coeff = F/v [N/(m/s)]
        m=0.5;  %   kg
        M=1;    %   kg
        l=0.6;  %   m
        g=9.81; %   m/s^2
        
        %srednja vrednost reference:
        r0=NaN;
        %PID regulator:
        Kpid=[100; 1;30;2];
    end
    
    methods
        function [Kpid]=getpidparam(self)
           Kpid=self.Kpid; 
        end
        function [r]=getrefsignal(self)
            dr=0.02;
            doff=0.10;

            %stevilo diskretnih vrednosti ene periode
            nzlepkov=3*2+3*2;
            np=floor(self.Tper/nzlepkov)/self.Ts*nzlepkov; 
            r1=zeros(np/nzlepkov,1);
            r21=repmat([r1-doff+dr; r1-doff-dr],[3,1])+self.r0;
            r22=repmat([r1+doff-dr; r1+doff+dr],[3,1])+self.r0;
            r2=[r21;r22];
            Np=floor(self.Tsim/self.Tper); %stevilo vseh period
            r=repmat(r2,[Np,1]);
        end
        function [u_ident]=getidentsignal(self)
            u_ident=randn(self.Tsim/self.Ts,1)*0.25;
        end
        function init(self)
            %sum
            self.noisestd=0e-3;
            %inicializacija stanj
            self.x.x=0;
            self.x.dx=0;
            self.x.fi=0;
            self.x.dfi=0;            
            %referencni signal se plete okrog r0:
            self.r0=0; 
        end
        function [y]=gety(self)
            yv=self.pendulumSim();
%             yv.x=yv.x+randn()*self.noisestd;
%             yv.fi=yv.fi+randn()*self.noisestd;
%             y=yv.x+self.l*sin(yv.fi)+randn()*self.noisestd;
            y=yv.fi+randn()*self.noisestd;
            self.y=y;
        end
        function setu(self,u)
%             u=self.saturate(u);
            self.u=u;
        end

        function self=pendulum_model(varargin)
            Ts=0.1;
            usat=struct('umin',-1,'umax',1,'ymin',-pi,'ymax',pi);
            normdata=struct('umean',0,'udelta',1,...
                            'ymean',0,'ydelta',0.5);
            if nargin>=1
                if ~(Ts==varargin{1})
                    error('Ts incorrect! Ts must be 0.5s');
                end
            end
            if nargin>=2
                usat=varargin{2};
            end
            if nargin==3
                normdata=varargin{3};
            end
            self=self@process_interface(Ts,usat,normdata);
        end
        
        function [y] = pendulumSim(self)
            f=@(t,x) self.odefun(t,x,self.u);
            x0=self.x;
            x0ode=[x0.x x0.dx x0.fi x0.dfi]';
            [t,x1ode]=ode23(f,[0 self.Ts],x0ode); %11 steps per each run is minimum
            x1.x   =x1ode(end,1);
            x1.dx  =x1ode(end,2);
            x1.fi  =x1ode(end,3);
            x1.dfi =x1ode(end,4);
            self.x=x1;
            %output:
            y.x=x1.x;
            y.fi=x1.fi;
        end
        function [ds]=odefun(self,t,s,u)
            ds=NaN(3,1);
            m=self.m;
            l=self.l;
            M=self.M;
            g=self.g;
            b=self.b;
            F=u;

            x =s(1);
            dx=s(2);
            fi =wrapTo2Pi(s(3)+pi)-pi;
            dfi =s(4);

            ds(1)=dx;
            ds(2)=(2*l*m*sin(fi)*dfi^2 + 4*F - 4*b*dx + ...
                3*g*m*cos(fi)*sin(fi))/(4*M + 4*m - 3*m*cos(fi)^2);
            ds(3)=dfi;
            ds(4)=-(3*(l*m*cos(fi)*sin(fi)*dfi^2 + 2*F*cos(fi) + 2*g*m*sin(fi) + ...
                2*M*g*sin(fi) - 2*b*dx*cos(fi)))/(l*(4*M + 4*m - 3*m*cos(fi)^2));
%               syms ddx ddfi
%               rez=solve(...
%               '(M+m)*ddx+m*l*ddfi*cos(fi)/2-m*l*dfi^2*sin(fi)/2==F-b*dx',...
%               '2*l*ddfi+3*ddx*cos(fi)+3*g*sin(fi)',...
%               ddx,ddfi)
        end
    end
end

