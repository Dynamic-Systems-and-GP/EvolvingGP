classdef EGP < matlab.mixin.Copyable

properties(GetAccess=public,SetAccess=public)

	BVi		=[];
	BVt		=[];
	BVtst	=[];
	
	hyp	;
	
	inf ;
	cov ;
	lik ;
	mean;
	size=0;
	logLikelihood=inf;
	post=struct();
	
	hypOptim;
	forgetting;
	reducing;
	signals;
end

methods
	function self=EGP(signals)
		self.signals=signals;
		
		self.inf = {'infExact'};
		self.cov = {'covSEard'};
		self.lik = {'likGauss'};
		self.mean= {'meanZero'};		
		
		self.hypOptim.enable=1;
		self.hypOptim.iter=10;
		self.forgetting.factor=1;
		self.forgetting.type='none';
		self.reducing.maxSize=100;
		self.reducing.type='window';
	end
	
	function resetPrior(self,default)
		if nargin<2 default=1; end
        if nargin==2 && isstruct(default)
           self.hyp=default;
           return;
        end
		%initialize default values:
		self.hyp=struct();
		D=self.signals.Nregressors;
		self.hyp.cov= ones(eval(feval(self.cov{:})),1)	.*default;
		self.hyp.mean=ones(eval(feval(self.mean{:})),1)	.*default;
		self.hyp.lik= ones(eval(feval(self.lik{:})),1)	.*default;
	end
	
	function inferPosterior(self)
		if isempty(self.BVi) return; end
		[self.post,nlZ]=feval(self.inf{:},self.hyp, self.mean, self.cov, self.lik{1}, self.BVi, self.BVt);
		self.logLikelihood=-nlZ;
	end
	
	function include(self,k)
        k=k(:);
        N=length(k);
		try
			x=self.signals.getRegressorVectors(k);
			t=self.signals.getRegressand(k);
		catch err
			if strcmp(err.identifier,'SignalsModel:InvalidSignalIndex')
                warning('regressors not avaliable.');
				return;
			else
				rethrow(err);
			end
		end
		self.BVi(end+1:end+N,:)=x;
		self.BVt(end+1:end+N,:)=t;
		self.BVtst(end+1:end+N,:)=k;
		self.size=size(self.BVi,1);
	end
	function [ymu,ys2]=predictAt(self,k)
        k=k(:);
        try
			xs=self.signals.getRegressorVectors(k);
		catch err
			if strcmp(err.identifier,'SignalsModel:InvalidSignalIndex')
				ymu=inf(size(k));
				ys2=inf(size(k));
			else
				rethrow(err);
			end
        end
        [ymu,ys2]=predict(self,xs);
    end
        
	function [ymu,ys2]=predict(self,xs)
        if self.size==0 ymu=Inf;ys2=Inf; return; end
		post=self.post;
		cov=self.cov;
		mean=self.mean;
		hyp=self.hyp;
		lik=self.lik;
		x=self.BVi;
		
		try
			%beware, the following code is partially modified from gp.m [Copyright (c) 
			% by Carl Edward Rasmussen and Hannes Nickisch, 2011-02-18]:
			  if isempty(lik),  lik = @likGauss; else                        % set default lik
				  if iscell(lik), lik = lik{1}; end                      % cell input is allowed
				  if ischar(lik), lik = str2func(lik); end        % convert into function handle
			  end
			  alpha = post.alpha; L = post.L; sW = post.sW;
			  if issparse(alpha)                  % handle things for sparse representations
				nz = alpha ~= 0;                                 % determine nonzero indices
				if issparse(L), L = full(L(nz,nz)); end      % convert L and sW if necessary
				if issparse(sW), sW = full(sW(nz)); end
			  else nz = true(size(alpha)); end                   % non-sparse representation
			  if numel(L)==0                      % in case L is not provided, we compute it
				K = feval(cov{:}, hyp.cov, x(nz,:));
				L = chol(eye(sum(nz))+sW*sW'.*K);
			  end
			  Ltril = all(all(tril(L,-1)==0));            % is L an upper triangular matrix?
			  ns = size(xs,1);                                       % number of data points
			  nperbatch = 1000;                       % number of data points per mini batch
			  nact = 0;                       % number of already processed test data points
			  ymu = zeros(ns,1); ys2 = ymu; fmu = ymu; fs2 = ymu;% lp = ymu;   % allocate mem
			  while nact<ns               % process minibatches of test cases to save memory
				id = (nact+1):min(nact+nperbatch,ns);               % data points to process
				kss = feval(cov{:}, hyp.cov, xs(id,:), 'diag');              % self-variance
			%     disp(kss);
				Ks  = feval(cov{:}, hyp.cov, x(nz,:), xs(id,:));         % cross-covariances
				ms = feval(mean{:}, hyp.mean, xs(id,:));
				fmu(id) = ms + Ks'*full(alpha(nz));                       % predictive means
				if Ltril           % L is triangular => use Cholesky parameters (alpha,sW,L)
				  V  = L'\(repmat(sW,1,length(id)).*Ks);
				  fs2(id) = kss - sum(V.*V,1)';                       % predictive variances
			%       disp(sum(V.*V,1)');
				else                % L is not triangular => use alternative parametrisation
				  fs2(id) = kss + sum(Ks.*(L*Ks),1)';                 % predictive variances
				end
				fs2(id) = max(fs2(id),0);   % remove numerical noise i.e. negative variances
				% if nargin<9
				  [~, ymu(id) ys2(id)] = lik(hyp.lik, [], fmu(id), fs2(id));
				% else
				  % [lp(id) ymu(id) ys2(id)] = lik(hyp.lik, ys(id), fmu(id), fs2(id));
				% end
				nact = id(end);          % set counter to index of last processed data point
			  end
		end
	end
	
	function optimizePrior(self)
		if isempty(self.BVi) return; end
		if self.hypOptim.enable==1
			[self.hyp, nlZ] = minimize(self.hyp, @gp, self.hypOptim.iter,...
				self.inf,self.mean, self.cov,  self.lik,self.BVi,self.BVt);
			self.inferPosterior();
		end
	end	
	
	function reduce(self,informationGain)
		if isempty(self.BVi) return; end
		exceededSize=self.size-self.reducing.maxSize;
		if exceededSize>0
			if nargin==1 informationGain=self.getInformationGain();	end
			informationGain=self.applyForgetting(informationGain);		
			timestamps=informationGain(end-exceededSize+1:end,2);
			id=NaN(exceededSize,1);
			for i=1:length(exceededSize)
				id(i)=find(self.BVtst==timestamps(i));
			end
			self.BVi(id,:)	=	[];
			self.BVt(id,:)	=	[];
			self.BVtst(id,:)=	[];		
			
			self.inferPosterior();
			%~ fprintf('reduced data timestamps: %s\n',num2str(timestamps(id)));				
			fprintf('                worsest information gain element is at k-%d step with value (%d)\n' ,max(self.BVtst)-informationGain(end,2),informationGain(end,1));
		end
		self.size=size(self.BVi,1);
	end
	
	function [informationGain]=getInformationGain(self)
		tst=self.BVtst;
		informationGain=NaN(self.size,2);
		informationGain(:,2)=tst(:);	
		
		switch self.reducing.type(1:3)
			case 'max'
				sig=1;
				reducingmethod=self.reducing.type(4:end);
			case 'min'
				sig=-1;
				reducingmethod=self.reducing.type(4:end);
			otherwise
				sig=1; % pretending to be "max", the preposition for maximum value
				reducingmethod=self.reducing.type(1:end);
		end
		
		
		
		switch reducingmethod
			case 'likelihood'
				  for i = 1 : self.size
						[~,nlZ] = feval(self.inf{:},self.hyp, self.mean, self.cov, self.lik{1}, self.BVi([1:i-1,i+1:end],:), self.BVt([1:i-1,i+1:end]));
						informationGain(i,1) = -nlZ*sig;
				  end
			case 'optimizedlikelihood'
				  for i = 1 : self.size
					[self.hyp, nlZ] = minimize(self.hyp, @gp, self.hypOptim.iter,...
						self.inf,self.mean, self.cov,  self.lik,self.BVi([1:i-1,i+1:end],:), self.BVt([1:i-1,i+1:end]));
						informationGain(i,1) = -nlZ(end)*sig;
				  end
			case 'euclid'
				% sorts the elements of active set by ascending euclid distance between them
				D=dist([self.BVi self.BVt]')*sig;
				D_id1=repmat(tst(1:self.size),1,self.size)';
				D_id2=D_id1';
				diagsizes=(self.size-1:-1:1);
				diagids=[0 cumsum(diagsizes)];
				D1=NaN(2*diagids(end),3);
				k=1;
				for i=1:length(diagsizes)
					D12((diagids(i)+1):diagids(i+1),:)=[diag(D,i),diag(D_id1,i),diag(D_id2,i)];
				end
				D1=[D12(:,[1,3]);D12(:,[1,2])];
				[~,D1i]=sort(D1(:,1),'descend');
				D1=D1(D1i,:);
				%~ if sig>0 occurrence='last'; else  occurrence='first'; end;
				occurrence='last';
				[informationGain(:,2),D1indices]=unique(D1(:,2),occurrence);
				informationGain(:,1)=D1(D1indices,1);
				%~ [min(D1(:,1)),min(informationGain(:,1))]
			case 'linearindependence'
			% sorts the elements of active set by ascending euclid distance between them
				[mu,se2]=self.predict(self.BVi);
				informationGain(:,1)=sqrt(se2);
			case 'windowing'
			% sorts the elements of active set by ascending euclid distance between them
				informationGain(:,1)=tst(:);
			otherwise
				error(['unknown reducing method: ' reducingmethod]);
		end
		informationGain=self.applyForgetting(informationGain);
	end
		
	function [informationGain] = applyForgetting(self,informationGain)
		%
		%
		% InformationGain is a (Nx2) vector of information gains with corresponding indeces.
		%
		% apply a forgetting (by sample age) term to an already computed 
		% information gain for each sample inside dataset. The forgetting 
		% factor is linear with correpsonding sample timestamps.
		
		%shift information gain to be equal or greater than 0
		informationGain(:,1)=informationGain(:,1)-min(informationGain(:,1));
		timestampnow=max(informationGain(:,2));
		switch self.forgetting.type
			case 'linear'
				s=informationGain(:,1)-self.forgetting.factor*(timestampnow-informationGain(:,2));
			case 'exponential'
				s=informationGain(:,1)*self.forgetting.factor^(timestampnow-informationGain(:,2));
			case 'none'
				s=informationGain(:,1);
		end
		[~,sid]=sort(s,'descend');
		informationGain=informationGain(sid,:);		
    end
    
end


end




	%~ egp.Ts=Ts; %->to dyn
    %~ egp.lambda=lambda; %-> MPC control
    %~ egp.eta=eta; %-> MPC control
    %~ egp.rho=rho; %-> MPC control
    %~ egp.rhoT=rhoT; %-> MPC control
    %~ egp.RegulatorEProfile=RegulatorEProfile; %-> MPC control
    %~ egp.RegulatorVProfile=RegulatorVProfile; %-> MPC control
    %~ egp.RegulatorUProfile=RegulatorUProfile; %-> MPC control
    %~ egp.uexp=uexp; %-> MPC control

    %~ egp.msteps=msteps;%% -> MPC control
    %~ egp.upd_steps=upd_steps; %% -> to dyn, specific for update condition
    
    %~ egp.udelta=process.udelta; %% lacks impl. for multi-input
    %~ egp.umean=process.umean; %% lacks impl. for multi-input
    %~ egp.ydelta=process.ydelta; %% lacks impl. for multi-output?
    %~ egp.ymean=process.ymean; %% lacks impl. for multi-output?
    %~ egp.BVt_id=[];%% -> change name to timestamps: egp.BVtst
