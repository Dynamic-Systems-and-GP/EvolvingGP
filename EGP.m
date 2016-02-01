classdef EGP < handle
    
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
            self.reducing.maxSize=NaN;
            self.reducing.type='windowing';
            self.reducing.indexLifes=NaN(signals.maxk,2);
        end
        
        function resetActiveSet(self,default)
            k=self.BVtst;
            self.BVi=[];
            self.BVt=[];
            self.BVtst=[];
            
            self.include(k);
            %self.inferPosterior();
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
            
            K = feval(self.cov{:},self.hyp.cov,self.BVi);
            L = chol(K + exp(2*self.hyp.lik)*eye(size(self.BVi,1)))';
            self.post.invQ = L'\(L\eye(size(self.BVi,1)));
            self.post.beta=L'\(L\self.BVt);
        end
        
        
        
        function include(self,k)
            k=k(:);
            
            [RVmask,RSmask]=self.signals.filterUndefinedDataAt(k,'warning');
            k=k(RVmask&RSmask);
            if isempty(k) return; end
            x=self.signals.getRegressorVectors(k);
            t=self.signals.getRegressand(k);
            N=length(k);
            
            self.BVi(end+1:end+N,:)=x;
            self.BVt(end+1:end+N,:)=t;
            self.BVtst(end+1:end+N,:)=k;
            self.size=size(self.BVi,1);
            
            for i=1:length(k)
                self.reducing.indexLifes(k(i),:)=[self.signals.time(k(i)),Inf];
            end
            
        end
        function [ymu,ys2]=predictAt(self,k,varargin)
            k=k(:);
            ymu=NaN(length(k),1);
            ys2=NaN(length(k),1);
            %         try
            RVmask=self.signals.filterUndefinedDataAt(k,'warning');
            newk=k(RVmask);
            %             ymu(~mask_k)=NaN;
            %             ys2(~mask_k)=NaN;
            [xs,vs]=self.signals.getRegressorVectors(newk);
            % 		catch err
            % 			if strcmp(err.identifier,'SignalsModel:InvalidSignalIndex')
            % 				ymu=inf(size(k));
            % 				ys2=inf(size(k));
            %                 return;
            % 			else
            % 				rethrow(err);
            % 			end
            %         end
            if any(any(vs))
                for i=1:length(newk)
                    maskids=find(RVmask);
                    [ymu(maskids(i)),ys2(maskids(i))]=propagate(self,xs(i,:),diag(vs(i,:)));
                end
            else
                [ymu(RVmask),ys2(RVmask)]=predict(self,xs);
            end
            
            if nargin==3
                if isempty(varargin{1})
                    
                elseif strcmpi(varargin{1},'normalized') || strcmpi(varargin{1},'n')
                    return;
                else
                    error(['unknown option: ' varargin{1}]);
                end
            end
            
            ymu=ymu*self.signals.std.(self.signals.O)+self.signals.mean.(self.signals.O);
            ys2=ys2*(self.signals.std.(self.signals.O)^2);
            
        end
        
        function [ymu,ys2]=predict(self,xs,varargin)
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
                    [~, ymu(id), ys2(id)] = lik(hyp.lik, [], fmu(id), fs2(id));
                    % else
                    % [lp(id) ymu(id) ys2(id)] = lik(hyp.lik, ys(id), fmu(id), fs2(id));
                    % end
                    nact = id(end);          % set counter to index of last processed data point
                end
            catch
                warning('prediction failed.');
            end
            
        end
        
        function optimizePrior(self)
            if isempty(self.BVi)
                return;
            end
            if self.hypOptim.enable==1
                [self.hyp, nlZ] = minimize(self.hyp, @gp, self.hypOptim.iter,...
                    self.inf,self.mean, self.cov,  self.lik,self.BVi,self.BVt);
                self.inferPosterior();
            else
                warning('Hyperparameter optimization is configured to be disabled! Nothing to be done. Leaving...');
            end
        end
        
        function reduce(self,informationGain)
            if isempty(self.BVi) return; end
            exceededSize=self.size-self.reducing.maxSize;
            if exceededSize>0
                if nargin==1
                    informationGain=self.getInformationGain();
                    informationGain=self.applyForgetting(informationGain);
                end
                timestamps=informationGain(end-exceededSize+1:end,2);
                id=NaN(exceededSize,1);
                for i=1:exceededSize
                    id(i)=find(self.BVtst==timestamps(i));
                end
                self.BVi(id,:)	=	[];
                self.BVt(id,:)	=	[];
                self.BVtst(id,:)=	[];
                
                self.inferPosterior();
                %~ fprintf('reduced data timestamps: %s\n',num2str(timestamps(id)));
                fprintf('                worsest information gain element is at k-%d step with value (%d)\n' ,max(self.BVtst)-informationGain(end,2),informationGain(end,1));
                
                for i=1:length(id)
                    self.reducing.indexLifes(id(i),2)=self.signals.time(id(i));
                end
                
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
                    % sorts the elements of active set by ascending linear
                    % independence in RKHS
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
                    s=informationGain(:,1)-self.forgetting.factor.*(timestampnow-informationGain(:,2));
                case 'exponential'
                    s=informationGain(:,1).*self.forgetting.factor.^(timestampnow-informationGain(:,2));
                case 'none'
                    s=informationGain(:,1);
            end
            [~,sid]=sort(s,'descend');
            informationGain=informationGain(sid,:);
        end
        function new = copy(self)
            % Copy super_prop
            new = feval(class(self),self.signals);
            
            % Copy all non-hidden properties.
            p = properties(self);
            for i = 1:length(p)
                new.(p{i}) = self.(p{i});
            end
            new.signals = self.signals.copy();
        end
        
        
        
        function [m, S2] = propagate(self, muX, SigX)
            % Computes the predictive mean and variance at teh test input with
            % Gaussian distribution for the squared exponential covariance function.
            % Input:
            % * hyp      ... optimized hyperparameters
            % * inf      ... the function specifying the inference method
            % * mean     ... the prior mean function
            % * cov      ... the specified covariance function, see help covFun for more info
            % * lik      ... the likelihood function
            % * invQ     ... the inverse of the data covariance matrix
            % * input    ... the input part of the training data,  NxD matrix
            % * target   ... the output part of the training data (ie. target), Nx1 vector
            % * muX      ... the 1 by D test input
            % * SigX     ... the covariance of the test input (OPTIONAL)
            %
            hyp=self.hyp;
            inf=self.inf;
            mean=self.mean;
            cov=self.cov;
            lik=self.lik;
            invQ=self.post.invQ;
            input=self.BVi;
            target=self.BVt;
            
            if iscell(cov) cov=str2func(cov{1}); end
            if iscell(inf) inf=str2func(inf{1}); end
            if iscell(lik) lik=str2func(lik{1}); end
            if iscell(mean) mean=str2func(mean{1}); end
            
            beta=self.post.beta;
            
            [n, D] = size(input); % the number of training cases and dimension of input space
            [nn, D] = size(muX);  % the number of test cases and dimension of input space
            
            % input validation
            
            [ is_valid, hyp, inf, mean, cov, lik, msg ] = validate( hyp, inf, mean, cov, lik, D);
            
            
            if ~isequal(cov,{@covSEard})
                error(strcat([fun_name,': function can only be called with the', ...
                    ' covariance function ''covSEard'' ']));
            end
            
            if ~isequal(lik,{@likGauss})
                error(strcat([fun_name,': function can only be called with the', ...
                    ' likelihood function ''likGauss'', where hyp.lik parameter is log(sn)']));
            end
            
            
            X=[-2*hyp.cov(1:end-1);2*hyp.cov(end);2*hyp.lik]; % adapt hyperparameters to local format
            expX = exp(X);        % exponentiate the hyperparameters once and for all
            
            %~ beta = invQ*target;
            
            % Covariance between training and test inputs ...
            
            a = zeros(n,nn);
            for d = 1:D
                a = a + (repmat(input(:,d),1,nn)-repmat(muX(:,d)',n,1)).^2*expX(d);
            end
            a = expX(D+1)*exp(-0.5*a);
            
            % Covariance between the test input and themselves
            b = expX(D+1);
            
            % Predictive mean  and variance (test input including noise variance)
            m = a'*beta;
            S2 = b - sum(a.*(invQ*a),1)'  + expX(D+2);
            
            if  nargin > 2 % nondeterministic test input
                
                L = sum(diag(SigX)~=0); % number of stochastic dimensions (non zero variances)
                if (L==0)
                    return;
                end
                % split regressor to deterministic quantities and
                % nondeterministic quantities
                rangeL = find(diag(SigX)~=0);
                rangeC = find(diag(SigX)==0);
%                 find(diag(SigX)~=0)
                
                SigX = SigX(rangeL,rangeL);
                muXL = muX(:,rangeL);
                muXC = muX(:,rangeC);
                inputL = input(:,rangeL);
                inputC = input(:,rangeC);
                
                invLL = diag(expX(rangeL));
                invLC = diag(expX(rangeC));
                invS = inv(SigX);
                invC = (invLL+invS);
                invSmuX = invS*muXL';
                t1 = muXL*invSmuX;
                c = inv(invC)*(invLL*inputL'+repmat(invSmuX,1,n));
                t2 = sum(inputL.*(inputL*invLL),2);
                t3 = sum(c.*(invC*c),1)';
                I = (1/sqrt(det(invLL*SigX+eye(L))))*exp(-0.5*(t1+t2-t3));
                CC = exp(-.5*(sum((inputC-repmat(muXC,n,1)).*((inputC-repmat(muXC,n,1))*invLC),2)));
                m = b.*(CC.*I)'*beta;
                
                invD = 2*invLL+invS;
                [kk1,kk2]=meshgrid(1:n);
                T1 = repmat(inputL,n,1)+inputL(reshape(kk1,1,n^2),:);
                invLT1 = T1*invLL;
                d = invD\(invLT1'+repmat(invSmuX,1,n^2));
                T3 = reshape(sum(d.*(invD*d),1),n,n);
                I2 = (1/sqrt(det(2*invLL*SigX+eye(L))))*exp(-0.5*(t1+repmat(t2,1,n)+repmat(t2',n,1)-T3));
                
                CCC = CC*CC';
                S2 = b - b^2*sum(sum((invQ-beta*beta').*(CCC.*I2))) - m^2 + expX(D+2); % including noise variance
                
            end
        end
        
        
    end
end