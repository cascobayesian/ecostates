library(VGAM)
library(gtools)


sim.data = function(N,alpha,counts)
{
	K = length(alpha);
	sims = array(NA,dim=c(N,K))

	for(i in 1:N)
	{
		draw = rdirichlet(1,alpha);	
		sims[i,] = rmultinom(1,counts,draw);
	}
	return(sims)
}


sample.v = function(data.set,alpha)
{
	z = 0;
	if (data.set>0)
	{	
		z = 1;
		if (data.set>1)
		{
			ms = seq(2,data.set,1);
			z = sum(c(z,rbinom(length(ms),1,alpha/(alpha+ms-1))));
		}
	}
	return(z)
}

sample.z = function(j,data.set,alpha)
{
	N = length(data.set);
	d = dim(as.matrix(data.set));

	#print(c(j,d))
	if (1 %in% d)
	{
		return(sample.v(data.set[j],alpha[j]))		
	}else
	{
		return(sum(unlist(lapply(data.set[,j],sample.v,alpha[j]))))
	}
}

sample.zs = function(data.set,alpha)
{
	K = length(alpha);
#	print(K)
	return(unlist(lapply(1:K,sample.z,data.set,alpha)))
}

gibbs.z = function(k,data.set,model)
{
	new.model = model;
	cs = which(new.model$c == k);
	new.model$Z[k,] = sample.zs(data.set[cs,],new.model$A[k,])	
	return(new.model)
}

gibbs.a = function(k,data.set,model)
{
	new.model=sample.lambda(k,model);
	new.model=sample.omega(k,data.set,new.model);
	return(new.model)
}

d.omega.post = function(w, Ns,model,N,k)
{
	val = dexp(w,0.1,log=TRUE) + sum(model$Z[k,])*log(w) + N*lgamma(w)-sum(lgamma(w+Ns))
	return(val)
}

sample.omega = function(k,data.set,model)
{
	cs = which(model$c==k)
	if (length(cs)>1)
	{
		Ns = rowSums(data.set[cs,]);
	}else
	{
		Ns = sum(data.set[cs,]);
	}
	ws = seq(0.0001,100*model$M,length.out=2000)
	
	y = unlist(lapply(ws,d.omega.post,Ns,model,length(cs),k));
	y = exp(y-max(y,na.rm=TRUE))
	
	new.model = model;
	new.model$O[k] = sample(ws,1,prob=y);

	new.model$A[k,] = new.model$O[k]*new.model$L[k,];
	return(new.model)
}

sample.lambda = function(k,model)
{
	new.model = model;
	new.model$L[k,] = rdirichlet(1,model$Z[k,]+1);
	return(new.model)
}

calc.llk.site = function(data.set,alpha)
{
	A = sum(alpha);
	N = sum(data.set);

	llk = lgamma(A) - lgamma(N+A) + sum(lgamma(data.set+alpha)) - sum(lgamma(alpha))

	return(llk)
}

calc.llk = function(data.set,model)
{
	llk = rep(0,model$N)
	
	for(i in 1:model$N)
	{
		llk[i] = calc.llk.site(data.set[i,],model$A[model$c[i],])
	}
	return(llk)
}

init.model = function(K, data.set, alpha=5)
{
	N = dim(data.set)[1]
	M = dim(data.set)[2];	
	c = rep(1:K,N)[1:N]
	
	A  = array(0,dim=c(K,M))
	Z  = A;
	L  = A;
	O  = rep(0,M);


	for(k in 1:K)
	{
		A[k,] = rep(1,M);
		Z[k,] = rep(5,M);
		L[k,] = A[k,]/sum(A[k,])
	}
	
	O  = rowSums(A)	
	model = list(c,A,Z,O,L,alpha,N,M)

	names(model) = c("c","A","Z","O","L","alpha","N","M")
	tables = calc.tables(model,data.set)		
	llk = calc.llk(data.set,model)
	
	model = list(c,A,Z,O,L,alpha,N,M,llk,tables[[1]],tables[[2]])

	names(model) = c("c","A","Z","O","L","alpha","N","M","llk","tables","counts")

	return(model)
}

calc.tables = function(model,data.set)
{
	tables = table(model$c)
	counts = array(0,dim=c(length(tables),model$M))

	clusts = as.numeric(names(tables))

	for(i in 1:length(tables))
	{
		these = which(model$c==clusts[i])
		counts[i,] = counts[i,] + colSums(data.set[these,])
	}
	return(list(tables,counts))
}


adjust.tables.2 = function(model,i,new.clust,new.a,m=3)
{

	new.model = model;
	current.clusters = as.numeric(names(model$tables))
	if (new.clust %in% current.clusters)
	{

		## then it's in here!
		this = which(current.clusters == model$c[i])
		new.model$tables[this] = new.model$tables[this]-1;
		this = which(current.clusters == new.clust)	
		new.model$tables[this] = new.model$tables[this]+1;	
		new.model$c[i] = new.clust;

	}else
	{

		this = which(current.clusters == model$c[i])
		new.model$tables[this] = new.model$tables[this]-1;
		new.model$tables = c(new.model$tables,1);
		names(new.model$tables) = c(names(model$tables),as.character(length(model$tables)+1))
		new.model$c[i] = length(new.model$tables);		
		new.model$A = rbind(model$A,new.a)
		new.model$O = c(model$O,sum(new.a));
		new.model$Z = rbind(model$Z,rep(1,model$M));
		new.model$L = rbind(model$L,new.a/sum(new.a));
	}

#	new.model$tables = table(new.model$c);
	zeros = which(new.model$tables==0)
	if(length(zeros)>0)
	{

		clusters  = as.numeric(names(new.model$tables))
		counter = 1
		for(i in 1:length(new.model$tables))
		{
			if (new.model$tables[i]>0)
			{
				these = which(new.model$c==clusters[i])
				new.model$c[these] = counter;
				counter = counter + 1;
			}
		}
		new.model$tables = new.model$tables[-zeros];
		new.model$A = new.model$A[-zeros,];
		new.model$Z = new.model$Z[-zeros,];
		new.model$L = new.model$L[-zeros,];
		new.model$O = new.model$O[-zeros]
	}

	names(new.model$tables) = c(1:length(new.model$tables))

	return(new.model)
}

draw.prior = function(M)
{
	return(rexp(M,1))
}

mh.draw.1 = function(i,model,data.set,m=3)
{
	current = model$c[i];
	tables = model$tables;
	counts = model$counts
	N      = model$N;
	M      = model$M
	alpha  = model$alpha;

	this = which(current == as.numeric(names(tables)))
	tables[this] = tables[this]-1;

	k.minus = length(tables[tables>0]);
	h = k.minus + m;

	tab.new = as.numeric(names(tables[tables>0]));
	
	new.A = array(0,dim=c(m,M))

	for(j in 1:m)
	{
		new.A[j,] = draw.prior(M)		
	}

	llks = rep(0,h);

	for(r in 1:h)
	{
		if (r <= k.minus)
		{
			llks[r] = calc.llk.site(data.set[i,],model$A[tab.new[r],]);			
		}else
		{
			llks[r] = calc.llk.site(data.set[i,],new.A[r-k.minus,]);
		}
	}

	probs = c(tables[tables>0],rep(alpha/m,m))/(N-1+alpha)*exp(llks-max(llks))
		
	new.c = sample(1:h,1,prob=probs);

	if (new.c>k.minus)
	{
		model = adjust.tables.2(model,i,new.c,rbind(model$A,new.A)[new.c,])
	}else
	{
		model = adjust.tables.2(model,i,new.c,NA)		
	}

#	print(c(length(model$tables),new.c,model$c[i]))
	model$llk[i] = calc.llk.site(data.set[i,],as.matrix(model$A[model$c[i],]))	
	return(model)

}


gibbs.draw = function(k,model,data.set)
{
	new.model = gibbs.z(k,data.set,model);
	new.model = gibbs.a(k,data.set,new.model);

	return(new.model)
}

mh.draw.2 = function(k,model,data.set)
{
	M = model$M;
	j = sample(1:M,1)
	A.new = model$A[k,]
	A.new[j] = draw.prior(1)	

	llk.new = rep(0,model$N)
	for(i in 1:model$N)
	{
		if (model$c[i] == k)
		{
			llk.new[i] = calc.llk.site(data.set[i,],A.new);
		}
	}	
	log.ratio = sum(llk.new)-sum(model$llk[model$c==k])
	log.u = log(runif(1))
	if(log.u < log.ratio)
	{
		model$A[k,] = A.new;
	}
	return(model);
}

mh.draw = function(model,data.set)
{
	N = model$N;
	for(i in 1:N)
	{
		model = mh.draw.1(i,model,data.set,m=3)
	}
	for(j in 1:length(model$tables))
	{
		model = gibbs.draw(j,model,data.set)
	}
	return(model)
}
