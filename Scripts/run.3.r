source("Scripts/model.3.r")
source("Scripts/dp.5.r")

args = commandArgs(trailingOnly=TRUE);

if (length(args)!=2)
{
	print("WTF?!?! You're done!" )
#	data.file = "Data/GOM.plankton/copepod.75000.tab"
#	iteration = 1
	q()
}else
{
	data.file = args[1] 
	iteration = as.numeric(args[2])
}

data.set = read.table(data.file,header=TRUE,sep="\t");
data.set = as.matrix(data.set);

num.iter = 3;
current = rep(list(),num.iter)
current[[1]] = init.model(15,data.set);

x = strsplit(data.file,"/")[[1]]
x = strsplit(x[3],"\\.")[[1]];

out = paste("Output/",x[1],".",x[2],".",iteration,".RData",sep="")
for(j in 2:num.iter)
{
	print(c(j,length(table(current[[j-1]]$c))))
	current[[j]] = mh.draw(current[[j-1]],data.set)
}

save(file = out,list=ls())

