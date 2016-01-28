% PROCESS -----------
%select a process to work with: (but it's not limited to the predefined process classes)
proc=bioreactor_icecream(0.5);
proc.noisestd=1e-3;
u=proc.getidentsignal();
proc.init();
N=length(u);

% select which regressors will be used for modelling:
%  create a struct where each property is a signal name with its regressors (the delayed signals)
inps.u=[1,3];
inps.y=[1,2];

% SIGNALS -----------
%Instantiate a 'signals' object called "es" and allocate him a memory by giving the length of our signal: N
es=SignalsModel(N);
%tell which regressors should we use:
%	The method setInputs creates signals inside the object es with same names as defined in "inps"
es.setInputs(inps);			
% tell which signal is the regressand (the target data)
es.setOutput('y');			

%since we already defined whole input signal u, we define it again inside the signals object:
%WARNING: CHANGED DEFINITION OF SIGNALS
% signals now contain two columns. one for mean, second for variances (default zero)...
es.u(:,1)=u(:);

%Allocate additional two signals to object 'es'. To store predicted model (mu, se2).
%We specifically want additional signals. They will be used for storing predictions of our model.
addprop(es,'yp_mu');
addprop(es,'yp_se2');

% EGP -----------
e=EGP(es);					%instantiate the EGP-model class

D=es.Nregressors;
hyp.cov= ones(eval(covSEard()),1);
hyp.lik=1;
hyp.mean=[];
e.resetPrior(hyp);				%initialize hyperparameters 
%e.resetPrior(100);				%initialize hyperparameters 
e.hypOptim.iter=5;			%set how much iterations should we do for opt.
e.reducing.maxSize=50;		%set the maximum size of active set

%the following property defines the method for calculation of information gain:
%	-the higher infromation gain for an excluded element, the lesser chance that
%   the element will be actually excluded from the active set.
%   -the information gains of elements follows the "min" or "max" property of the selected methods:
%		1.) min/max Euclid distance of i-the element to any other  (method: '(min|max)euclid')
%		2.) min/max marginal likelihood by excluding i-th element from active set  (method: '(min|max)likelihood')
%		3.) min/max marginal optimized likelihood w.r.t. the hyperparameters. It is similar to (2.) 
% 			but runs an optimization for each element   (method: '(min|max)optimizedlikelihood')
%		4.) the youthness of each element (method:windowing)
%
%		Hint: if the method is named without the min/max prefix, the default is max value of the method result.
e.reducing.type='maxeuclid';%what should the information gain of a single (active set) element equal to.

for k=1:N
	es.k=k; %optional: 		tell the 'signals' object which is the current time-step
	
	es.y(k)=proc.gety();	%simulate the process and get the current output
    if k>50 && k<60 | k>400 && k<470
        es.y(k)=NaN;					% let y(k) cannot be measured for 51<k<59 and 401<k<470
    end	
	[es.yp_mu(k), es.yp_se2(k)]=e.predict(k-1);	%predict the same output with the model
	
	% if abs(error)>tolerance then update the model with new regressor 
	%							vector and regressand at time-step "k"
	if abs(es.yp_mu(k)-es.y(k))>0.005	
		e.include(k);				% increase the active set by one
		e.optimizePrior();			% optimize hyps
		e.inferPosterior();			% calculate all vectors for fast prediction, e.g., the covariance inverse
		e.reduce();					% reduce the active set. The method verifies the size of active set by itself
	end
	proc.setu(es.u(k,1));				% apply the current input signal value to the process
end

ax(1)=subplot(2,1,1);
plot(1:N,es.y(:,1),'b',1:N-1,es.yp_mu(2:end),'g');
ax(2)=subplot(2,1,2);
plot(es.u(:,1));
linkaxes(ax,'x');

proc.delete();
