classdef bioreaktor_model_recycle < process_interface
    %RAKETA_PROCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess=public,SetAccess=public)
        Tsim=500;
        Tper=70;
        %realnocasovne simulacije ne bomo strogo uporabljali s timestamp
        timestamp=NaN;
        %vhodi, izhodi, stanje sistema
        u=NaN; % this is F -- outflow rate
        y=NaN; % this is X, the biomass
        yr=NaN;
        S=NaN;
        X=NaN;
        
        %brezdimenzijski parametri modela [Nelson et at. 2008]:
		P=struct('tau',NaN,'R',NaN,'k_d',NaN,'X_0',NaN,'alpha',NaN,'m_s',NaN,'beta',NaN,'gamma',NaN);
        %~ P=struct([]);
        
        %srednja vrednost reference:
        r0=0.4;
        
        tau=0;
        quantize=0;
        %PID regulator:
        Kpid=[-14; -5; 0;.5];
    end
    
    methods
        function [Kpid]=getpidparam(self)
           Kpid=self.Kpid; 
        end
        
    function [r]=getrefsignal(self,reftype)
        switch lower(reftype)
        %type twotrains
            case 'twotrain'
                dr=0.05;
                doff=0.15;

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
                Nlvls=round(self.Tsim/(self.Ts*10));
                Ntimes=floor(self.Tsim/(Nlvls*self.Ts));

                r=randn(1,Nlvls);%/5+0.4;
                r=self.ymean+r*self.ydelta*1/5;
                r=reshape(repmat(r,[Ntimes,1]),[Nlvls*Ntimes,1]);
            
            case 'constant'
                r=ones(round(self.Tsim/self.Ts),1)*self.r0;
        end
    end 
        function [u_ident]=getidentsignal(self)
            Tch=10;
            Nch=200;
            ur=rand(1,Nch);%/5+0.4;
            ur=ur*(self.umax-self.umin)+self.umin;
            u_ident=reshape(repmat(ur,[Tch/self.Ts,1]),[Nch*Tch/self.Ts,1]);
        end
        function init(self,varargin)
            %~ self.P.tau		=1;			% dimensionless residence time (input signal)
			self.P.R		=0.9;		% effective recycle parameter
			self.P.k_d		=0.1;		% dimensionless decay rate (from death rate)
			self.P.X_0		=0;			% dimensionless concentration of microorganisms in the influent
			self.P.alpha	=1;			% dimensionless yield coefficient
			self.P.m_s		=0.04;		% dimensionless maintenance energy
			self.P.gamma	=1;%sem si izmislil
			self.P.beta		=1;%sem si izmislil
        
            %inicializacija stanj
            self.S=0.8;
            self.X=0.2;
            self.y=self.X;
            self.u=self.umean;
            self.yr=self.X;
        end
        function [y]=gety(self)
            self.bioreactorSim();
            y = self.X+randn*self.ydelta*self.noisestd;
            self.y=y;
        end
        function [yr]=getyr(self)
            yr = self.X;
        end
        function u=setu(self,u)
            u=self.saturate(u);
            self.u=u;
        end
		
		function bioreactorSim(self)
		    f=@(t,s) self.odefun(t,s,self.u);
            s0ode=[self.S self.X]';
            [t,s1ode]=ode23(f,[0 self.Ts],s0ode); %11 steps per each run is minimum
            self.S   =s1ode(end,1);
            self.X   =s1ode(end,2);
		end
		
        function [ds]=odefun(self,t,s,u)
            ds=NaN(2,1);
			
			R		=self.P.R		;
			k_d		=self.P.k_d		;
			X_0		=self.P.X_0		;
			m_s		=self.P.m_s		;
			alpha	=self.P.alpha	;
			beta	=self.P.beta	;
			gamma	=self.P.gamma	;
            
            % dimensionless residence time tau' = V/(F*mu_max), F is the
            % real input (u) and let's assume V/mu_max = 1 [1/g]
            tau=1/u; % assume t
            S=s(1);
			X=s(2);
			
            ds(1)=(1/tau)*(1-S) - (1/alpha)*(X*S)/(X+S) - m_s*X;				% dS/dt
            ds(2)=beta*(1/tau)*(X_0-X) + gamma*R/tau*X + (X*S)/(X+S) - k_d*X;	% dX/dt
        end
        
        function [u]=findEquilibriumInput(self)
            u=self.umean;
        end  
        
        function self=bioreaktor_model_recycle(varargin)
            Ts=0.5; % in HOURS = 1800s
%             usat=struct('umin',0.05,'umax',1,'ymin',0.03,'ymax',1.34);
            usat=struct('umin',0.2,'umax',10,'ymin',0.03,'ymax',1.34);
            normdata=struct('umean',(usat.umax+usat.umin)/2,'udelta',abs(usat.umax-usat.umin),...
                            'ymean',0.688,'ydelta',0.3);
%             normdata=struct('umean',0.35,'udelta',0.35,...
%                             'ymean',0.688,'ydelta',0.3);
            if nargin>=2
                usat=varargin{2};
            end
            if nargin==3
                normdata=varargin{3};
            end
            self=self@process_interface(Ts,usat,normdata);
        end
    end
end

