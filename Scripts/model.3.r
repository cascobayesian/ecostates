initialize.model = function(data.set,K,var=10)
{
	require(mvtnorm)
	
	num.sites = dim(data.set)[1];
	num.count = dim(data.set)[2];

	Z = rep(1:K,num.sites)[1:num.sites];

	beta  = array(0,dim=c(K,num.count))
	P     = array(NA,dim=c(K,num.count));

	for(i in 1:K)
	{
		P[i,]     = exp(beta[i,])/sum(exp(beta[i,]))
	}

	sigma= diag(rep(var,num.count));
	mu   = rep(0,num.count);

	llk = 0;
	lpr = 0;
	
	model = list(num.sites,num.count,K,Z, P,beta,sigma,mu,llk,lpr,0,0)
	names(model) = c("M","D","K","Z","P","beta","sigma","mu","llk","lpr","llk.vec","lq")

	model$llk.mat = calc.llk.mat(data.set,model)	
	model$llk = calc.llk(model)
	model$lpr = calc.lpr(model);
	return(model)
}


calc.llk.mat = function(data.set,model)
{
	llk.mat= array(NA,dim=c(model$M,model$K))
	for(j in 1:model$M)
	{
		for(k in 1:model$K)
		{
			llk.mat[j,] = dmultinom(data.mat[j,],size=sum(data.mat[j,]),prob=model$P[k,],log=TRUE)
		}
	}	
	return(llk.mat)
}

calc.llk = function(model)
{
	llk = 0;
	for(j in 1:model$M)
	{
		llk = llk+model$llk.mat[j,model$Z[j]];
	}
	return(llk)
}

calc.lpr = function(model)
{
	lpr = 0
	for(k in 1:model$K)
	{
		lpr = lpr + dmvnorm(model$P[k,],mean=model$mu,sigma=model$sigma,log=TRUE)
	}
	return(lpr)
}

propose.z = function(model,data.set)
{
	js = sample(1:model$M,floor(model$M/4),replace=FALSE)
	new.model =model;
	new.model$lq=0;
	for(j in js)
	{

		llks = model$llk.mat[j,]
		probs = exp(llks-max(llks))
	
		cur.k = model$Z[j];
		new.k = sample(1:model$K,1,prob=probs)

		new.model$Z[j] = new.k;
		new.model$lq  = new.model$lq+log(probs[cur.k])-log(probs[new.k]);
	}
	new.model$llk = calc.llk(new.model)

	return(new.model)
}

propose.b = function(model,data.set)
{
	new.k = sample(1:model$K,1);
	new.d = sample(1:(model$D-1),1);
	new.model = model;
	new.model$beta[new.k,new.d] = rnorm(1,mean=0,sd=1);

	new.model$P[new.k,] = exp(new.model$beta[new.k,])/sum(exp(new.model$beta[new.k,]))

	for(j in 1:model$M)
	{
		new.model$llk.mat[j,new.k] = dmultinom(data.mat[j,],size=sum(data.mat[j,]),prob=new.model$P[new.k,],log=TRUE)
	}	
	new.model$llk = calc.llk(new.model);
	new.model$lpr = calc.lpr(new.model);	

	new.model$lq  = dnorm(model$beta[new.k,new.d],mean=0,sd=1,log=TRUE) - dnorm(new.model$beta[new.k,new.d],mean=0,sd=1,log=TRUE);

	return(new.model)
}


draw = function()
{
	opts = c('b','z')
	d = sample(opts,1,prob=c(0.8,0.2))
	return(d);	
}

propose.model=function(data.set,model)
{
	d = draw();
	if (d=='b')
	{
		new.model=propose.b(model,data.set)
	}else
	{
		new.model=propose.z(model,data.set)
	}	
	return(new.model)
}

mcmc = function(data.set,K,num.iter=30000,thin=1000)
{
	current = test = initialize.model(data.set,K)	
	chain = rep(list(),num.iter/thin);

	for(i in 1:num.iter)
	{
		proposed = propose.model(data.set,current)
	
		log.p = proposed$lpr - current$lpr;
		log.l = proposed$llk - current$llk;
		log.q = proposed$lq;

		log.a = log.p+log.l+log.q;
		if (log(runif(1))< log.a)
		{
			current = proposed;
		}
		if (i %% thin == 0)
		{
			chain[[i/thin]]=current;
			print(c(i/thin,current$llk));
		}
	}
	return(chain)
}
