classdef linunstable_model < process_interface
    %RAKETA_PROCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess=public,SetAccess=public)
        Tsim=301;
        Tper=150;
        tau=0.5;
        %realnocasovne simulacije ne bomo strogo uporabljali s timestamp
        timestamp=NaN;
        %parametri iz modela:
        u=NaN;
        y=NaN;
        x=NaN;
        %srednja vrednost reference:
        r0=0;
        %PID regulator:
        Kpid=[6; 0.02; 10;0.4];
        quantize=NaN;
        
        A=[];
        B=[];
        C=[];
        D=[];
        
        mcont=[];
        mdisc=[];
    end
    
    methods
        function [Kpid]=getpidparam(self)
           Kpid=self.Kpid; 
        end
        function [r]=getrefsignal(self,reftype)
            switch lower(reftype)
            %type twotrains
                case 'twotrain'
                    dr=0.03;
                    doff=0.07;

                    %stevilo diskretnih vrednosti ene periode
                    nzlepkov=3*2+3*2;
                    np=floor(self.Tper/nzlepkov)/self.Ts*nzlepkov; 
                    r1=zeros(np/nzlepkov,1);
                    r21=repmat([r1-doff+dr; r1-doff-dr],[3,1])+self.r0;
                    r22=repmat([r1+doff-dr; r1+doff+dr],[3,1])+self.r0;
                    r2=[r21;r22];
                    Np=floor(self.Tsim/self.Tper); %stevilo vseh period
                    r=repmat(r2,[Np,1]);

            %type random
                case 'random'
                    Nlvls=round(self.Tsim/(self.Ts*10*5));
                    Ntimes=floor(self.Tsim/(Nlvls*self.Ts));

                    r=randn(1,Nlvls);%/5+0.4;
                    r=self.ymean+r*self.ydelta*1/5;
                    r=reshape(repmat(r,[Ntimes,1]),[Nlvls*Ntimes,1]);

                case 'constant'
                    r=ones(round(self.Tsim/self.Ts),1)*self.r0;
            end
        end
            
        function [r]=getrefsignal0(self)
            dr=0.005;
            doff=0;%0.01; 
        
            %stevilo diskretnih vrednosti ene periode
            nzlepkov=12;
            np=floor(self.Tper/nzlepkov)/self.Ts*nzlepkov; 
            r1=zeros(np/nzlepkov,1);
            r21=repmat([r1-doff+dr; r1-doff-dr],[3,1])+self.r0;
            r22=repmat([r1+doff-dr; r1+doff+dr],[3,1])+self.r0;
            r2=[r21;r22];
            Np=floor(self.Tsim/self.Tper); %stevilo vseh period
            r=repmat(r2,[Np,1]);
            r=[rand(np/nzlepkov/2,1)*doff+self.r0     ;   r];
        end
        function [u_ident]=getidentsignal(self)
            Tch=4;
            Nch=100;
            ur=rand(1,Nch);%/5+0.4;
            ur=ur*(self.umax-self.umin)-self.umin;
            u_ident=reshape(repmat(ur,[Tch/self.Ts,1]),[Nch*Tch/self.Ts,1]);
        end
        function init(self,varargin)
            
            p=-1*2*pi;
            ceta=0.2;
            w=-0.5;
            
            self.mcont=tf([1],[1,2*ceta*w,w^2])*tf([1],[1 -p]);
            self.mdisc=ss(c2d(self.mcont,self.Ts));
            self.A=self.mdisc.A;
            self.B=self.mdisc.B;
            self.C=self.mdisc.C;
            self.D=self.mdisc.D;
            
            

            %inicializacija stanj
            self.x=zeros(size(self.C'));
            self.y=self.C*self.x;
            self.u=0;
            %referencni signal se plete okrog r0:
            self.r0=0;
        end
        function [y]=gety(self)
            self.x = self.A*self.x+self.B*self.u;
            y = self.C*self.x+self.D*self.u+randn*self.noisestd*self.ydelta;
            self.y=y;
        end
        function u=setu(self,u)
            u=self.saturate(u);
            self.u=u;
        end

       
        function self=linunstable_model(varargin)
            Ts=0.1;
            usat=struct('umin',-0.5,'umax',0.5,'ymin',-0.25,'ymax',0.25);
            normdata=struct('umean',0,'udelta',1,...
                            'ymean',0,'ydelta',0.5);
            if nargin>=1
                Ts=varargin{1};
            end
            if nargin>=2
                usat=varargin{2};
            end
            if nargin==3
                normdata=varargin{3};
            end
            self=self@process_interface(Ts,usat,normdata);
        end
        
        function [xout] = narendra_p(self, u,x, varargin)       
            xout=x/(1+x^2)+u^3+randn*self.noisestd; 
        end
%         function [dx, y] = narendrali_m(t, x, u, p, varargin)
        %NARENDRALI_M  A discrete-time Narendra-Li benchmark system.
        
        
        
        %   Copyright 2005-2006 The MathWorks, Inc.
        %   $Revision: 1.1.8.1 $ $Date: 2006/11/17 13:28:56 $

        % Output equation.
%         y = x(1)/(1+p(4)*sin(x(2)))+x(2)/(1+p(5)*sin(x(1)));
% 
%         % State equations.
%         dx = [(x(1)/(1+x(1)^2)+p(1))*sin(x(2));              ... % State 1.
%               x(2)*cos(x(2))+x(1)*exp(-(x(1)^2+x(2)^2)/p(2)) ... % State 2.
%                  + u(1)^3/(1+u(1)^2+p(3)*cos(x(1)+x(2)))     ...
%              ];
%         
    end
end







