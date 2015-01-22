classdef SignalsModel < dynamicprops & matlab.mixin.Copyable

properties(GetAccess=public,SetAccess=public)	
	time;
	k;
	maxk;
	maxlag;
	Nregressors;
end

properties
	I={};
	Ilags={};
	O='';
	Ninputs;
	Noutputs;
end

methods
	function self=SignalsModel(N)
		self.maxk=N;
		self.time=NaN(N,1);
		self.k=1;
		%~ Nn=length(signalnames)
		%~ for i=1:Nn
			%~ addprop(self,signalnames{i});
			%~ self.set(signalnames{i})=NaN(N,1);
		%~ end
	end
	
	function [x]=getRegressorVectors(self,k)
		if nargin==1 k=self.k; end;
		assert(min(k)>self.maxlag,'SignalsModel:InvalidSignalIndex','Too small time-step. The regressor vector requires lagged values at/before time-step k=0');
		k=k(:);
		n=length(k);
		x=NaN(n,self.Nregressors);
		j=1;
		for i=1:self.Ninputs
			nr=numel(self.Ilags{i});
			x(:,j:j+nr-1)=self.(self.I{i})(bsxfun(@minus,k,self.Ilags{i}));
			j=j+nr;
        end
        if any(isnan(x))
            error('SignalsModel:InvalidSignalIndex','Invalid regressor values.');
        end
		
	end
	function [t]=getRegressand(self,k)
		if nargin==1 k=self.k; end;
		k=k(:);
		n=length(k);
		t=self.(self.O)(k);
        
        if any(isnan(t))
            error('SignalsModel:InvalidSignalIndex','Invalid regressor values.');
        end
	end	
	function setInputs(self,inputs)
		%	setInputs(self,(struct) inputs)
		%
		% set the input signals for modelling. The argument "inputs" is a struct containing attributes named as the signals names. The value of each (attribute) name contains a vector that desribes which delayed signal values (l1,..,ln) from current time-step k are used for building a regressor vector.
		%
		% example:
		%	inputs.u=[1,3,5];
		%	inputs.y=[1,2];
		%	%in this case the feature vector (of current time-step k) for prediction is: x(k)=[u(k-1) u(k-3) u(k-5) y(k-1) y(k-2)];
		%	"obj".setInputs(inputs);
		signalnames=fieldnames(inputs);
		signaldelays=struct2cell(inputs);
		Nin=length(signalnames);
		for i=1:Nin
			addprop(self,signalnames{i});
			self.I{i}=signalnames{i};
			self.Ilags{i}=signaldelays{i};
			self.(signalnames{i})=NaN(self.maxk,1);
		end
		self.Ninputs=Nin;
		self.maxlag=max([self.Ilags{:}]);
		self.Nregressors=numel([self.Ilags{:}]);
	end
	function setOutput(self,output)
		%	setOutput(self,(str) output)
		%
		% set the output (target) signal for modelling
		if ~isprop(self,output)
			addprop(self,output);
			self.(output)=NaN(self.maxk,1);
		end
		self.O=output;
		self.Noutputs=1;
    end
    
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
end
