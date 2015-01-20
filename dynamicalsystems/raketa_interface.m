classdef raketa_interface < process_interface
    %PROCESS_INTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    
properties(GetAccess=public,SetAccess=public)
        %cas
        Tsim=1500;
        Tper=280;%Tper=130;
        
        tau=2;
%         Tsim=210;
%         Tper=200;

        timestamp=NaN;
        %sum
%         noisestd=0;%1e-2;
        
        %ref signal
        r0=0.285;
        
        %parametri iz modela:
        H=0.71;         %m
%         y0=0.3; %range -5 to 1.146
        u=NaN;
        y=NaN;
        x0=struct('z',NaN,'dz',NaN,'p',NaN);
        x=struct('z',NaN,'dz',NaN,'p',NaN);
        %PID regulator:
        Kpid=NaN(1,4);
end

methods
    function [r]=getrefsignal(self,reftype)
        switch lower(reftype)
        %type twotrains
            case 'twotrain'
                dr=0.02;
                doff=0.03;  
                Nvagonov=5;

                %stevilo diskretnih vrednosti ene periode
                nzlepkov=Nvagonov*2+Nvagonov*2;
                np=floor(self.Tper/nzlepkov)/self.Ts*nzlepkov; 
                r1=zeros(np/nzlepkov,1);
                r21=repmat([r1-doff+dr; r1-doff-dr],[Nvagonov,1])+self.r0;
                r22=repmat([r1+doff-dr; r1+doff+dr],[Nvagonov,1])+self.r0;
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
    function self=raketa_interface(varargin)
        Ts=0.1;
%         usat=struct('umin',0,'umax',5,'ymin',0.046,'ymax',0.48);
%         normdata=struct('umean',2.25,'udelta',1,...
%                         'ymean',0.285,'ydelta',0.185);
        usat=struct('umin',1,'umax',3,'ymin',0.046,'ymax',0.48);
        normdata=struct('umean',2.25,'udelta',1,...
                        'ymean',0.285,'ydelta',0.185);  
                    
                    %original udelta=2.5
                    %umean je 2.25 ko je raketa polna zraka
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
        
        %imena in enote signalov
        self.ulabel='u [V]';
        self.ylabel='h [m]';
    end
    function init(self)
        %Not implemented
    end
    function [y]=gety(self)
        %Not implemented
    end
    function setu(self,u)
        %Not implemented
    end
    function [usat]=saturate(self,u)
        if u>self.umax u=self.umax
        elseif u<self.umin u=self.umin; end
        usat=u;
    end
    function delete(self)
        %Not implemented
    end
end

end