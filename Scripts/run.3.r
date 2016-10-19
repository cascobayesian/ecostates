#!~/anaconda/bin/Rscript
# NOTE: CHANGE THE ABOVE TO YOUR RSCRIPT PATH
# =============================================================================
# =============================================================================
#                              ecostates v1.2
# =============================================================================

# =============================================================================
# Infinite dimensional Dirichlet-multinomial mixtures for ecological count data
# =============================================================================
# GitHub: https://github.com/jacobian1980/ecostates
# Publication: http://biorxiv.org/content/early/2016/04/20/045468
# Authors: John D. O'Brien, Nicholas Record, Peter Countway
# doi: http://dx.doi.org/10.1101/045468
# =============================================================================



# =============================================================================
# Edit Log:
# =============================================================================
# Josh L. Espinoza (github.com/jolespin)
# ______________________________________
#	08.22.2016
#		i . Added label functionality (Now takes in data with labels.  First column must be row labels)
#	09.30.2016
# 		ii . Added seed for reproducibility
# 		iii . Produces output for all iterations in 3 DataFrames (For cluster assignment, log-likelihood, and probability vectors)
#		iv . Added shebang for terminal execution
#	10.04.2016
#		v . Changed argument length error
#		vi . Fix sourcing when not in current directory
#		vii . Cleaned stdout
#   10.07.2016
#       viii . Added dimensions to stderr
#       ix . Added average log-likelihood for each iteration
# =============================================================================
# Pending
# =============================================================================
#	i . Add parameter for default starting cluster number
#	ii . Add option for stopping algorithm upon convergence
#	iii . Exit properly




# =============================================================================
# Utility Functions
# =============================================================================
library(base)
getExecutablePath <- function() {
        args = commandArgs(trailingOnly = FALSE)
        executable = "--file="
        file_path = grep(executable, args)
        if (length(match) > 0) {
                return(normalizePath(sub(executable, "", args[file_path])))
        }
}
executable_path = dirname(getExecutablePath())

# =============================================================================
# Import
# =============================================================================
cat("\n=============================================================================\necostates_v1.2\n=============================================================================\n")
cat("====Sourcing scripts...\n\n")
source(paste(executable_path, "model.3.r", sep="/"))
source(paste(executable_path, "dp.6.r", sep="/"))


# =============================================================================
# Set random seed
# =============================================================================
set.seed(0)


# =============================================================================
# Arguments
# =============================================================================
args = commandArgs(trailingOnly=TRUE);
if (length(args)!=2)
{
	print("Not enough arguments.  Did you add OTU counts table and number of iterations?" )
	q()
}else
{
	data.file = args[1]
	iteration = as.numeric(args[2])
}


# =============================================================================
# Load datasets
# =============================================================================
data.set = read.table(data.file,header=TRUE,sep="\t");

row_labels = data.set[[1]]
data.set = data.set[,-1]
data.set = as.matrix(data.set);

num.iter = iteration;
current = rep(list(),num.iter)
current[[1]] = init.model(15,data.set);

x = tail(strsplit(data.file,"/")[[1]], 1)

out = paste(x,".",iteration,".RData",sep="")


# =============================================================================
# Parameters
# =============================================================================
cat("====Parameters...\n================\n")
cat(paste("Input File:", data.file, "\n", sep="\t"))
cat(paste("Output Directory:", getwd(), "\n", sep="\t"))
cat(paste("Iterations:", num.iter, "\n", sep="\t"))
cat(paste("Dimensions:", nrow(data.set), ncol(data.set), "\n", sep="\t"))
cat(paste("Date:", Sys.Date(), "\n\n", sep="\t"))

# =============================================================================
# Run Algorithm
# =============================================================================
cat("====Run algorithm...\n=================\n")
cat(paste("iteration","num_clusters", "mu(llk)", "\n", sep="\t"))
for(j in 2:num.iter)
{
	cat(paste(j, length(table(current[[j-1]]$c)), mean(current[[j-1]]$llk), "\n", sep="\t"))
	current[[j]] = mh.draw(current[[j-1]],data.set)
}
cat("\n====Saving...\n")


# =============================================================================
# Organize Data
# =============================================================================
# Cluster | Log-likelihood Trace
cluster_trace = list(); likelihood_trace = list()
for (j in 1:num.iter){
	cluster_trace[[j]] = current[[j]]$c
	likelihood_trace[[j]] = current[[j]]$llk
	}
# Probability Trace
prob_data = list()
for (j in 1:num.iter){
	num_clusters = length(unique(current[[j]]$c))
	DF_tmp = as.data.frame(current[[j]]$L)
	# Add clusters
	DF_tmp = cbind(clusters = 1:num_clusters, DF_tmp, row.names = NULL)
	# Add iteration
	DF_tmp = cbind(iteration = rep(c(j), num_clusters), DF_tmp, row.names = NULL)
	colnames(DF_tmp) = c("iteration","cluster", colnames(data.set))
	prob_data[[j]] = DF_tmp
	}

# DataFrames
DF_c = as.data.frame(t(data.frame(matrix(unlist(cluster_trace), nrow=num.iter, byrow=TRUE), stringsAsFactors=FALSE)))
DF_c = cbind(obsv_ids = row_labels, DF_c, row.names = NULL)
colnames(DF_c) = c("obsv_id", 1:num.iter)

DF_llk = as.data.frame(t(data.frame(matrix(unlist(likelihood_trace), nrow=num.iter, byrow=TRUE), stringsAsFactors=FALSE)))
DF_llk = cbind(obsv_ids = row_labels, DF_llk, row.names = NULL)
colnames(DF_llk) = c("obsv_id", 1:num.iter)

DF_L = do.call("rbind", prob_data)


# =============================================================================
# Save Data
# =============================================================================
# Save DataFrames
write.table(DF_c, paste(x, ".", num.iter,  ".", "clusters", ".tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
write.table(DF_L, paste(x, ".", num.iter, ".", "probs", ".tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
write.table(DF_llk, paste(x, ".", num.iter, ".", "llk", ".tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")

# Save RData
save(file = out,list=ls())

# =============================================================================
# Complete
# =============================================================================
cat("====Complete...\n===============\n")
q()
