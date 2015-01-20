classdef raketa_proces < raketa_interface
    %RAKETA_PROCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess=public,SetAccess=public)
    end
    
    methods
        function self=raketa_proces(varargin)
            Ts=0.01;
            if nargin>=1
                Ts=varargin{1};
            end
            self=self@raketa_interface(Ts);
            %PID regulator:
%             self.Kpid=[-0.9/1.1; -.1; -1.5*3;0.08]/2; 
%             self.Kpid=[-0.9; -.1; -1.2;0.1]; 

            %Discrete PID regulator1: MAE=0.8933
%             self.Kpid=[-0.9; -.1; -1.5;0.1];     
            %Discrete PID regulator1: MAE=0.7567
%             self.Kpid=[-0.9; -.1; -1.2;0.1]; 
            %Discrete PID regulator2: MAE=0.4833
%             self.Kpid=[-0.8; -.1; -0.8;0.1];
%             Discrete PID regulator2: MAE=0.4833
            self.Kpid=[[-0.8; -.1; -0.8]*4.386 ;0.1];
        end
        function [u_ident]=getidentsignal(self)
            ur0=randn(self.Tsim/self.Ts,1)*0.25;
            u_ident=dlsim([1/4],[2,-1.9],ur0);
        end
        function init(self)
            % Real process drivers:
            if libisloaded('libdaqraketa')==1
                unloadlibrary('libdaqraketa');
            end
            loadlibrary('libdaqraketa','libdaqraketa.h');
            if calllib('libdaqraketa','isinit')==0
               calllib('libdaqraketa','init'); 
            end
        end
        function y=gety(self)
            y=calllib('libdaqraketa','readvalue');
            y= ((y-0.3)/2.5)*0.57;
        end
        
        function setu(self,u)
            u=self.saturate(u);
            calllib('libdaqraketa','sendvalue',u);
        end
        function [Kpid]=getpidparam(self)
           Kpid=self.Kpid;
        end
        function delete(self)
            calllib('libdaqraketa','sendvalue',0.0);
            if calllib('libdaqraketa','isinit')==1
                calllib('libdaqraketa','deinit');
            end
        end
    end
end

