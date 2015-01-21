classdef bioreaktor_model_icecream < process_interface
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
		P=struct('V',NaN,'Y',NaN,'B',NaN,'k_s',NaN,'k_d',NaN,'mu_max',NaN,'X_0',NaN);
        %~ P=struct([]);
        
        %srednja vrednost reference:
        r0=0.6;
        
        tau=0;
        quantize=0;
        %PID regulator:
        Kpid=[-5; -1; 0;.5];
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
            Nch=100;
            ur=rand(1,Nch);%/5+0.4;
            %~ ur=ur*(self.umax-self.umin)+self.umin;
            
            Du=(self.umax-self.umin)/3;
            Au=(self.umax+self.umin)/2;
            ur=ur*Du+Au;
            
            u_ident=reshape(repmat(ur,[Tch/self.Ts,1]),[Nch*Tch/self.Ts,1]);
        end
        function init(self,varargin)
            self.P.V		=5;			% bioreactor volume [l]
			self.P.Y        =0.2116;	% yield coefficient [ (g_VSS l^-1) / (g_COD l^-1) ] = [-]
			self.P.B		=0.4818;	% kinetic parameter (Contois type - K_s equivalent for Monod type) [g_COD / g_VSS]
			self.P.K_s		=0.4028;	% kinetic parameter (Monod type, not used here...) [g_COD]
            self.P.k_d		=0.0131;	% death rate [day^-1]
            self.P.mu_max	=0.9297;	% max population growth rate (Contois type) [day^-1]
			self.P.X_0		=0;			% biomass in the influent [g_VSS l^-1]
			self.P.S_0		=4.94  ;	% substrate in the influent [g_COD l^-1]
			
        
            %inicializacija stanj
            self.S=0.5075;
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
		
        function [ds_dt]=odefun(self,t,s,u)
            ds_dt=NaN(2,1);
			
			V		=self.P.V		;
			Y       =self.P.Y       ;
			B		=self.P.B		;
			K_s		=self.P.K_s		;
			k_d		=self.P.k_d		;
			mu_max	=self.P.mu_max	;
			X_0		=self.P.X_0		;
            S_0		=self.P.S_0		;
            
            Q=u;
            S=s(1);
			X=s(2);
			
			mu		 = mu_max* (S)/(B*X+S);
            ds_dt(1) = Q/V * (S_0 - S)  -  (1/Y)*X*mu ;			% dS/dt
            ds_dt(2) = Q/V * (X_0 - X)  +  mu*X - k_d*X;		% dX/dt
        end
        
        function [u]=findEquilibriumInput(self)
            u=self.umean;
        end  
        
        function self=bioreaktor_model_icecream(varargin)
            Ts=0.5; % in HOURS = 1800s
            usat=struct('umin',0,'umax',7,'ymin',0.3,'ymax',1.5);
            normdata=struct('umean',3.5,'udelta',3.5,...
                            'ymean',1,'ydelta',1.5);
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

%hot test hotline:
%clear classes;Ts=0.3;p=bioreaktor_model_icecream(Ts);p.init();p.gety();u=p.getidentsignal();N=length(u);y=NaN(N,1);for i=1:N y(i)=p.gety(); p.setu(u(i)); end; ax(1)=subplot(2,1,1); t=1:N; plot(t,y);ax(2)=subplot(2,1,2); plot(t,u); linkaxes(ax,'x');
