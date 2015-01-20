classdef bioreaktor_model < process_interface
    %RAKETA_PROCES Summary of this class goes here
    %   Detailed explanation goes here
    
    %Lastnosti objekta Bioreaktorja
    properties(GetAccess=public,SetAccess=public)
        Tsim=600;
        Tper=70;
        %realnocasovne simulacije ne bomo strogo uporabljali s timestamp
        timestamp=NaN;
        %parametri iz modela:
        u=NaN;
        y=NaN;
        x1=NaN;
        x2=NaN;
        %srednja vrednost reference:
        r0=0.07;
        tau=0;
        quantize=0;
        %PID regulator:
        Kpid=[-10; 0; 0;.5];
    end
    
    methods
        function [Kpid]=getpidparam(self)
            Kpid=self.Kpid;
        end
        
        %Referenèni signali
        function [r]=getrefsignal(self,reftype)
            switch lower(reftype)
                %type twotrains
                case 'twotrain'
                    dr=0.005;
                    doff=0.014;
                    
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
        
        %Vh. signal za identifikacijo
        function [u_ident]=getidentsignal(self)
            Tch=4;
            Nch=100;
            ur=rand(1,Nch);%/5+0.4;
            ur=ur*(self.umax-self.umin)-self.umin;
            u_ident=reshape(repmat(ur,[Tch/self.Ts,1]),[Nch*Tch/self.Ts,1]);
        end
        
        %Inicializacija sistema (bioreaktorja)
        function init(self,varargin)
            %inicializacija stanj
            self.x1=0.08;
            self.x2=0.02;
            self.y=self.x1;
            self.u=self.umean;
            %referencni signal se plete okrog r0:
            self.r0=0.066; %0.07 pri diplomi
            self.noisestd=0.01;
        end
        
        %Simuliraj in vrni izhodni signal
        function [y]=gety(self)
            [self.x1,self.x2] = self.bioeq(self.u,self.x1,self.x2);
            y = self.x1+randn*self.ydelta*self.noisestd;
            self.y=y;
        end
        
        %Vrni izhodni signal brez šuma
        function [yr]=getyr(self)
            yr = self.x1;
        end
        
        %Nastavi vhodni (regulirni) signal
        function u=setu(self,u)
            u=self.saturate(u);
            self.u=u;
        end
        
        %Matematièni model bioreaktorja
        function [x1n,x2n] = bioeq(self,uk,x1k,x2k)
            x1n = x1k + 0.5*x1k*x2k/(x1k+x2k) - 0.5*uk*x1k;
            x2n = x2k - 0.5*x1k*x2k/(x1k+x2k) - 0.5*uk*x2k + 0.05*uk;
        end
        
        %Konstruktor
        function self=bioreaktor_model(varargin)
            Ts=0.5;
            usat=struct('umin',0,'umax',0.7,'ymin',0.03,'ymax',0.1);
            normdata=struct('umean',0.35,'udelta',0.35,...
                'ymean',0.066,'ydelta',(0.1-0.03)/2);
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
    end
end

