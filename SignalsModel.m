classdef SignalsModel < dynamicprops

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
    mean;
    std;
    normalize;
end

methods
	function self=SignalsModel(N)
		self.maxk=N;
		self.time=NaN(N,1);
		self.k=1;
        self.normalize=1;
        self.mean=struct;
        self.std=struct;
		%~ Nn=length(signalnames)
		%~ for i=1:Nn
			%~ addprop(self,signalnames{i});
			%~ self.set(signalnames{i})=NaN(N,1);
		%~ end
	end
	
	function [x,v]=getRegressorVectors(self,k,varargin)
		if nargin==1 k=self.k; end;
        opt=0;
        if (nargin==3 && strcmpi(varargin{1},'force')) opt=1; end 
        if isempty(k) error('SignalsModel:getRegressorVectors','Empty time-step vector k'); return; end;
        
		if opt==0
           assert(min(k)>self.maxlag,'SignalsModel:InvalidSignalIndex','Too small time-instant. The regressor vector requires lagged values time-instant k=>0');
        else
           k(k<=self.maxlag)=[]; 
        end
        
		k=k(:);
		n=length(k);
		x=NaN(n,self.Nregressors);
		v=NaN(n,self.Nregressors);
		j=1;
        inputs_without_regressors=cellfun(@isempty,self.Ilags);
        for i=1:self.Ninputs
            if inputs_without_regressors(i)==1   continue;       end
            nr=numel(self.Ilags{i});
            x(:,j:j+nr-1)=self.(self.I{i})(bsxfun(@minus,k,self.Ilags{i}),1);
            v(:,j:j+nr-1)=self.(self.I{i})(bsxfun(@minus,k,self.Ilags{i}),2);
            if self.normalize==1
                x(:,j:j+nr-1)=(x(:,j:j+nr-1)-self.mean.(self.I{i}))./self.std.(self.I{i});
                v(:,j:j+nr-1)= v(:,j:j+nr-1)./(self.std.(self.I{i}).^2);
            end
            j=j+nr;
        end
        if any(any(isnan(x))) && opt==0
            error('SignalsModel:InvalidSignalIndex','Invalid regressor values.');
        end
    end
    
	function [t]=getRegressand(self,k)
		if nargin==1 k=self.k; end;
		k=k(:);
		n=length(k);
		t=self.(self.O)(k,1);
        if self.normalize==1
            t=(t-self.mean.(self.O))./self.std.(self.O);
        end
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
			self.(signalnames{i})=[NaN(self.maxk,1), zeros(self.maxk,1)]; % first column is mean, second is variance
            self.mean.(signalnames{i})=0;
            self.std.(signalnames{i}) =1;
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
			self.(output)=NaN(self.maxk,2);
            self.mean.(output)=0;
            self.std.(output) =1;
		end
		self.O=output;
		self.Noutputs=1;
    end
    
    function [regressorVector_mask,regressand_mask]=filterUndefinedDataAt(self,k,varargin)
        regressorVector_mask=[];
        regressand_mask=[];
       if nargin==3
           warntype=varargin{1}; % warning, error or quiet
       else
           warntype='warning';
       end
       k=k(:);
       
        shortk=(k<=self.maxlag);
        if any(shortk)
           self.(warntype)('SignalsModel:InvalidSignalIndex','Too small time-step. The regressor vector requires lagged values at/before time-step k=0');
        end
        
        k1=k(~shortk);
        if isempty(k1)
           k2=[]; %=k
           return;
        end
        
		n=length(k1);
		x=NaN(n,self.Nregressors);
		j=1;
        inputs_without_regressors=cellfun(@isempty,self.Ilags);
        for i=1:self.Ninputs
            if (inputs_without_regressors(i)==1)   continue; end
            nr=numel(self.Ilags{i});
            x(:,j:j+nr-1)=self.(self.I{i})(bsxfun(@minus,k1,self.Ilags{i}),1);
            j=j+nr;
        end
        o=self.(self.O)(k1,1);
        xsum=sum(x,2);
        nanx=isnan(xsum);
        nano=isnan(o);
        if any(nanx)>0 || any(nano)>0
            self.(warntype)('SignalsModel:InvalidSignalIndex','Invalid regressor values.');
        end
        k_mask=true(size(k));
        k_mask(shortk)=false;
        regressorVector_mask=k_mask;
        regressand_mask=k_mask;
        regressorVector_mask(~shortk)=~nanx;
        regressand_mask(~shortk)=~nano;
    end
    
    function updateStatMoments(self)
       %select only signal segments without NaN values
       for i=1:self.Ninputs
           notNaNs=~isnan(self.(self.I{i}),1);
           if sum(notNaNs)>1
             self.mean.(self.I{i})=mean(self.(self.I{i})(notNaNs,1));
             self.std.(self.I{i})=std(self.(self.I{i})(notNaNs,1));
           end
       end
    end
    
    function new = copy(self)
         % Copy super_prop
         new = feval(class(self),self.maxk);
         % Copy sub_prop1 in subclass
         % Assignment can introduce side effects
         new.setInputs(cell2struct(self.Ilags,self.I,2));
         new.setOutput(self.O);
         
         for i=1:length(self.I)
            signame=self.I{i};
            new.(signame)=self.(signame);
         end
         new.time=self.time;
		 new.k=self.k;
         new.mean=self.mean;
         new.std=self.std;
         new.normalize=self.normalize;
    end
    
    function quiet(self,id,msg)
        %empty function: is a quiet function, it isn't?
    end
    function warning(self,id,msg)
        warning(id,msg);
    end
    function error(self,id,msg)
        error(id,msg);
    end
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
end
