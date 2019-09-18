import argparse
import random
from random  import choice
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


class effectSizes:

        @classmethod
        def effectSizeList(cls,a,d):
                effectList={}
                for i in range(0, len(a)):
                        effectList[i]={"2":a[i],"1":d,"0":-a[i]}

                return effectList


class BaseHaplotypes:
	def __init__(self, ne):
		self.ne=ne
		
		if self.ne<0:
			raise ValueError("Number of diploid individuals should be larger than 0")

	@classmethod
	def make_basepop(self, popsize, counts):
		basepop=[]
		
		for l in counts:
			basepop.append(random.sample([1 for i in range(0,l)] + [0 for i in range(l,popsize)], popsize))

		diploids=make_diploids(basepop,popsize/2)
		return diploids




class make_diploids:
	def __init__(self,haplotypes,ne):
		self.haplotypes=haplotypes
		self.ne=ne

	def convert_to_diploids(self):
		diploid_basepop=[]
		for locus in range(0, len(self.haplotypes)):
			
			x=self.haplotypes[locus][::2]
			y=self.haplotypes[locus][1::2]

			diploid_basepop.append(list(zip(x,y)))

		diploid_basepop=[list(list(zip(*diploid_basepop))[i]) for i in range(0,int(self.ne))]
		return(diploid_basepop)


class individuals:
	def __init__(self,genotype,phenotype, fitness):
		self.genotype=genotype
		self.phenotype=phenotype
		self.fitness=fitness

	def geno(self):
		return self.genotype
	def pheno(self):
		return self.phenotype
	def fit(self):
		return self.fitness

	def gametes(self):
		gams=[]
		for i in self.genotype:
			gams.append(random.choice(i))
		return(gams)

class PopInfo:
	def __init__(self, population,effects,qff):
		self.population=population
		self.effects=effects
		self.qff=qff

	def popFunc(self):
		mypop=[]
		for ind in self.population.convert_to_diploids():
			phenotype=getPhenotype(ind,self.effects).estimation()
			fitness=getFitness(phenotype,self.qff).estimation()
			ind=individuals(ind, phenotype,fitness)
			mypop.append(ind)
		return(mypop)



class getPhenotype:
	def __init__(self,genotype,effect):
		self.genotype=genotype
		self.effects=effects

		if len(self.genotype)!=len(self.effects):
			raise ValueError("Inconsistent number of selected loci and effect sizes")				
	
	def estimation(self):
		ph=0
		for locus in range(0, len(self.genotype)):
			ph+=(self.effects[locus][str(sum(self.genotype[locus]))])
		return(ph)
			



class getFitness:
	def __init__(self, phenotype,qff):
		self.phenotype=phenotype
		self.minFit=qff[0]
		self.maxFit=qff[1]
		self.mean=qff[2]
		self.sd=qff[3]
	
	def estimation(self):
		scl=norm.pdf(self.mean, loc=self.mean, scale=self.sd)
		fitness=norm.pdf(self.phenotype, loc=self.mean, scale=self.sd)
		diff=self.maxFit-self.minFit
		sc_fitness=(fitness*diff/scl) +self.minFit
		return sc_fitness
	




class Selection:
	def __init__(self, pop,effects,qff):
		self.pop=pop
		self.popsize=len(self.pop)
		self.effects=effects
		self.qff=qff
		
	def fitnessTuple(self):
		ft=[]
		evoPop=[]
		ftsum=0
		
		for i in self.pop:
			ftsum+=i.fit()
			ft.append((ftsum,i))				
		for i in range(0,self.popsize):
			r1=random.uniform(0,ftsum)
			r2=random.uniform(0,ftsum)			
			ind1,f1=binarySearch(ft,r1)
			ind2,f2=binarySearch(ft,r2)		
			g1=ind1.gametes()
			g2=ind2.gametes()
			genotype=list(zip(g1,g2))
			phenotype=getPhenotype(genotype,self.effects).estimation()
			fitness=getFitness(phenotype,self.qff).estimation()
			ind=individuals(genotype, phenotype,fitness)
			evoPop.append(ind)

		return(evoPop)



def binarySearch(fitTuple,rand,lo=0, hi=None):
	
	if hi==None:
		hi=len(fitTuple)	

	while lo<hi:
		mid=int((lo+hi)/2)
		midval=fitTuple[mid][0]
		if midval<rand:			
			lo=mid+1
		elif(midval>rand):
			hi=mid
		else:
			return (fitTuple[mid+1][1],mid+1)
	return (fitTuple[lo][1],lo)





def getPopPheno(population):
	meanPheno=0
	phenoList=[]
	for pop in population:
		meanPheno=0
		for ind in pop[2]:
			meanPheno+=ind.pheno()
		avgPheno=meanPheno/len(pop[2])
		phenoList.append(avgPheno)


	return(phenoList)

def plotPhenotypes(populations,qff):
	mean=qff[2]
	box_plot_data=[]
	phenos=getPopPheno(populations)	
	
	if args.gen<=20:
		for i in range(0,args.gen+1):
#		box_plot_data.append(phenos[((args.gen+1)*i):(args.gen+((args.gen+1)*i))+1])
			box_plot_data.append(phenos[i::args.gen+1])
	
#		print(len(box_plot_data), len(list(map(str,range(0,args.gen)))))
		plt.boxplot(box_plot_data,labels=list(map(str,range(0,args.gen+1))))
		plt.axhline(y=mean, linewidth=1,linestyle='--', color = 'black')

		plt.show()
	
	elif args.gen>20 and args.gen <=100:
		for i in range(0,args.gen+1,5):
			box_plot_data.append(phenos[i::args.gen+1])
		plt.boxplot(box_plot_data,labels=list(map(str,range(0,args.gen+1,5))) )
		plt.axhline(y=mean, linewidth=1,linestyle='--', color = 'black')
		plt.xticks(rotation=90)
		plt.show()

	else:
		for i in range(0,args.gen+1,10):
			box_plot_data.append(phenos[i::args.gen+1])
		plt.boxplot(box_plot_data,labels=list(map(str,range(0,args.gen+1,10))) )
		plt.axhline(y=mean, linewidth=1,linestyle='--', color = 'black')
		plt.xticks(rotation=90)
		plt.show()




parser=argparse.ArgumentParser(description= """
            Description
            -----------
            Python script that simulates the evolutionary trajectory of a quantitative trait under stabilizing selection
	    (for diploids).

            Authors
            -----------
            Vlachos Christos""",formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument("-Ne",type=int, required=True,dest="ne",default=None, help="Population Size")
parser.add_argument("-a",type=str, required=True,dest="a",default=None, help="List of effect sizes for every QTL")
parser.add_argument("-p",type=str, required=True,dest="p",default=None, help="Initial Allele Frequency of every QTL")
parser.add_argument("-d",type=float, required=True,dest="d",default=None, help="Dominance effect")
parser.add_argument("-qff",type=str, required=True,dest="qff",default=None, help="Quantitative fitness function in the form: minFitness:maxFitness:mean:sd")
parser.add_argument("--replicates",type=int, required=True,dest="repl",default=None, help="Number of replicates")
parser.add_argument("--generations",type=int, required=True,dest="gen",default=None, help="Number of generations")
#parser.add_argument("--drift-distribution", required=False,dest="distr",default=False, help="If True returns the distribution of allele frequencies in the last generation (s must be 0)")
args = parser.parse_args()

popsize=args.ne*2
pcounts= list(map(int, [i*popsize for i in list(map(float,args.p.split(',')))]))
a=list(map(float, args.a.split(',')))
d=args.d
qff=list(map(float, list(args.qff.split(':'))))
effects=effectSizes.effectSizeList(a,d)
allpops=[]
print(effects)

for rep in range(1, int(args.repl)+1):
	basepop=BaseHaplotypes.make_basepop(popsize,pcounts)
	population=PopInfo(basepop,effects,qff).popFunc()
	allpops.append((rep,0,population))

	for g in range(1,int(args.gen)+1):
		print("Simulating generation {0} of replicate {1}".format(g,rep))
		population=(Selection(population,effects,qff).fitnessTuple())		
		allpops.append((rep,g,population))
	print("\n")

#print(allpops)
plotPhenotypes(allpops,qff)


#pop.selected_pop
