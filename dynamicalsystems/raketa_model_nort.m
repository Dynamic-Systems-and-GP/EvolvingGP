classdef raketa_model_nort < raketa_interface
    %RAKETA_PROCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess=public,SetAccess=public)
        
        quantize=1;
        yr=NaN;
    end
    
    methods
        function self=raketa_model_nort(varargin)
            Ts=0.1;
            if nargin>=1
                Ts=varargin{1};
            end
            self=self@raketa_interface(Ts);
            addpath('/home/martin/Dropbox/diploma/raketa_model');
            
%             PID regulator:
%             self.Kpid=[-.9*1.2,  -.1,    -1.5*2, 0.05]; 
%             self.Kpid=[-0.9, -.05, -1.5,0.1];
%               self.Kpid=[-3.5, -.1, -7,0.3];
            %Discrete PID regulator:
%             self.Kpid=[-3, -.1, -5,0.2];  
              a=2.2;
              b=2;
%               self.Kpid=[[-1.9, -0.07/a, -2.6]*a,0.15];            
%               self.Kpid=[[-1.2, -0.6/a, -2.2]*a,0.08];          
%                 self.Kpid=[[-0.8, -0.6/a, -2.4*b]*a,0.14*b];  
                self.Kpid=[-2.9, -0.5, -6.,0.18]; 
%             self.Kpid=[[-1; -.05; -0.5]*4.386;0.1];
        end
        
        function [u_ident]=getidentsignal(self)
            ur0=randn(self.Tsim/self.Ts,1)*0.25;
            u_ident=dlsim([1/4],[2,-1.9],ur0);
        end
        function init(self,varargin)
            %zacetek stetja casa:
            self.timestamp=clock();
            
            %parametri iz modela:
            self.H=0.57;
            H=0.57;
            m = 0.09075;    %kg
            S = 3.8e-4;     %m^2
            l = 0.145;      %m
            Vr = 95.85e-6;   %m^3
            g = 9.81;       %m/s^2
            ro = 1000;  %kg/m^3
            %labilna lega, zacetni pogoji - fiksiramo y0 -> z0:
            
            %inicializacija stanj
            if nargin>1
                z0=self.H-varargin{1};
            else
                z0=0.05;
            end
            dz0=0;
            
            %(dz0 mora biti nic, da ne bo dinamike)
            % y0=H-0.055;
            h0=1/(ro*S)*(ro*Vr-m); %delta visine
            f=@(p)abs(nivoVode(z0,p)-h0)*1e3; %iz funkcije poiscemo se tlak (numericno)
            p0 = fminbnd(f,0.5e5,2e5);

            self.x0=struct('z',z0,'dz',dz0,'p',p0);
            self.x=struct('z',z0,'dz',dz0,'p',p0);
            self.u=0;
            self.y=varargin{1};
            self.yr=varargin{1};
            
            if isprop(self,'noisestd')==0
                warning('setting noisestd to 1e-3.')
                self.noisestd=1e-3;
            end
            if self.noisestd>0
                warning('Noise is artificially added to process/model.')
            end
        end
        function [y,ysim]=gety(self)
            [ysim,x]=raketaSim(self.x,self.u,self.Ts);
            self.yr=ysim;
            y=ysim   +   randn()*self.noisestd*self.ydelta;
            if self.quantize==1
                y=quant(y,0.57/(64-7)); %kvantiziramo
            end
            % (64 - 7) kvantov, ker imamo 64 LEDic in zadnjih sedem ne bo raketa
            % nikoli dosegla.
            
            self.x=x;
            self.y=y;
        end
        function [yr]=getyr(self)
            yr=self.yr;
        end
        function [Kpid]=getpidparam(self)
           Kpid=self.Kpid;
        end
        function setu(self,u)
            u=self.saturate(u);
            self.u=u;
        end
        function [u]=findEquilibriumInput(self)
            ftmp=@(u) (raketaSim(self.x,u,self.Ts)-self.yr)^2*1e8;
            u=fminunc(ftmp,self.umean);
        end
%         function delete(self)
%         % Nothing to be done
%         end
    end
end

