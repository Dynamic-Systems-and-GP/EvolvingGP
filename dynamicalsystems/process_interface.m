classdef process_interface < handle
    %PROCESS_INTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    
properties(GetAccess=public,SetAccess=public)
    Ts=NaN;
    
    umin=NaN;
    umax=NaN;
    ymin=NaN;
    ymax=NaN;
    
    umean=NaN;
    udelta=NaN;
    ymean=NaN;
    ydelta=NaN;
    
    %imena in enote signalov
    tlabel='t [s]';
    ulabel='u';
    ylabel='y';
    
    noisestd=NaN;
end

methods
    function self=process_interface(Ts,sat,normdata)
        self.Ts=Ts;
        self.umin=sat.umin;
        self.umax=sat.umax;
        self.ymin=sat.ymin;
        self.ymax=sat.ymax;
        self.umean=normdata.umean;
        self.udelta=normdata.udelta;
        self.ymean=normdata.ymean;
        self.ydelta=normdata.ydelta;
    end
    function new = copy(this)
        % Instantiate new object of the same class.
        new = feval(class(this));

        % Copy all non-hidden properties.
        p = properties(this);
        for i = 1:length(p)
            new.(p{i}) = this.(p{i});
        end
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
            cond=[u<self.umin (u>=self.umin && u<=self.umax) u>self.umax];
            usat=cond*[self.umin u self.umax]';
            if ~cond(2)
                warning(['Input signal is over limits, saturating u=' num2str(usat) ]);
            end
    end
    function delete(self)
        %Not implemented
    end
end

end