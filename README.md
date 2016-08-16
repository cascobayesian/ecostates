# ecostates
Infinite dimensional Dirichlet-multinomial mixtures for ecological count data

This is a brief introduction to use the infinite dimensional Dirichlet-multinomial mixture model to understand ecological count data. The method is explained in detail [here](http://dx.doi.org/10.1101/045468). Essentially, if you have a data set of samples, where each sample is a set of counts across a set of taxa the method will infer the number of underlying components (distributions of taxa) and how the various samples are composed. It will also give you a posterior distribution on the number of components, via a Dirichlet process. 

This model is implemented in a set of R scripts in the above folders. The first step is to download these into an appropriate folder on your computer. Put whatever data you would like to use, in the data. To run it directly, you'll be working in the directory above that (so that if you folder is *foo* then your data is in *foo/Data*, but your commands will be in from *foo*.) Your data file will need to have the following features:
* the file should have a header but no row names;
* the file should be tab separated; and
* each row is a sample; each column is a taxa. 

To run the scripts, in your working directory type: 

``Rscript Scripts/run.3.r Data/file.txt num.iter &`` 

where ``file.txt`` is the name of your data file and ``num.iter`` is the number of iterations desired. Be warned: like all Dirichlet process mixture models, the inference can become sticky (i.e. converges poorly) as the number of samples grows and so larger number of iterations need to be used with more samples. The number of taxa alters this a bit, but the number of samples is the most important determinant. The ampersand will force the computer to treat the script as a stand-alone process. 

Once complete, the output is saved into the working directory as ``file.name.num.iter.RData``, where ``.RData`` is an R output format. This gives you the state of the Markov chain for each iteration. You then need to use additional scripts to process the chain and make sense of the data. The chain is an object named ``current`` and the model components are all located within this framework. The output can then be processed using additional scripts.
