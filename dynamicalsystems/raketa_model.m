classdef raketa_model < raketa_interface
    %RAKETA_PROCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess=public,SetAccess=public)
    end
    methods
        function self=raketa_model(varargin)
            Ts=0.1;
            if nargin>=1
                Ts=varargin{1};
            end
            self=self@raketa_interface(Ts);
            addpath('/home/martin/Dropbox/diploma/raketa_model');
            
            %PID regulator:
%             self.Kpid=[-.9*1.2,  -.1,    -1.5*2, 0.05]; 
            self.Kpid=[[-.9, -.1, -1.5]*4.386,0.1];
%             self.Kpid=[-0.8999   -0.0983   -1.5034    0.0621];
        end
        function [u_ident]=getidentsignal(self)
            ur0=randn(self.Tsim/self.Ts,1)*0.25;
            u_ident=dlsim([1/4],[2,-1.9],ur0);
        end
        function init(self)
            %zacetek stetja casa:
            self.timestamp=clock();
            
            %inicializacija stanj
            z0=0.05;
            dz0=0;
            
            %parametri iz modela:
            m = 0.09075;    %kg
            S = 3.8e-4;     %m^2
            l = 0.145;      %m
            Vr = 95.85e-6*2;   %m^3
            g = 9.81;       %m/s^2
            ro = 1000;  %kg/m^3
            %labilna lega, zacetni pogoji - fiksiramo y0 -> z0:
            %(dz0 mora biti nic, da ne bo dinamike)
            % y0=H-0.055;
            h0=1/(ro*S)*(ro*Vr-m); %delta visine
            f=@(p)abs(nivoVode(z0,p)-h0); %iz funkcije poiscemo se tlak (numericno)
            p0 = fminbnd(f,0.5e5,2e5);

            self.x0=struct('z',z0,'dz',dz0,'p',p0);
            self.x=struct('z',z0,'dz',dz0,'p',p0);
            self.u=0;
            self.y=0.57;
        end
        function [y]=gety(self)
            Telapsed=etime(clock,self.timestamp);
            [y,x]=raketaSim(self.x,self.u,Telapsed,self.noisestd);
            self.x=x;
            self.y=y;
            self.timestamp=clock();
        end
        function [Kpid]=getpidparam(self)
           Kpid=self.Kpid;
        end
        function setu(self,u)
            u=self.saturate(u);
            self.u=u;
            Telapsed=etime(clock,self.timestamp);
            [y,x]=raketaSim(self.x,self.u,Telapsed,self.noisestd);
            self.x=x;
            self.y=y;
            self.timestamp=clock();
        end
%         function delete(self)
%         % Nothing to be done
%         end
    end
end

