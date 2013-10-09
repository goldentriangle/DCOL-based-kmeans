""" 
Author: Kai Wang
Email: kai.wang.magic@gmail.com
"""
from scipy import *
from pylab import *
from sys import exit
from math import fabs
from random import random
import sys

def Distance(R1, R2):
	if len(R1)!=len(R2):
		print 'vector dimension disagree'
		exit()
	return sqrt( sum ( [(R1[i]-R2[i])**2 for i in range(len(R1))])  ) 

def TotalDistance(city, R):
    dist=0
    for i in range(len(city)-1):
        dist += Distance(R[city[i]],R[city[i+1]])
    dist += Distance(R[city[-1]],R[city[0]])
    return dist
    
def reverse(city, n):
    nct = len(city)
    nn = (1+ ((n[1]-n[0]) % nct))/2 # half the lenght of the segment to be reversed
    # the segment is reversed in the following way n[0]<->n[1], n[0]+1<->n[1]-1, n[0]+2<->n[1]-2,...
    # Start at the ends of the segment and swap pairs of cities, moving towards the center.
    for j in range(nn):
        k = (n[0]+j) % nct
        l = (n[1]-j) % nct
        (city[k],city[l]) = (city[l],city[k])  # swap
    
def transpt(city, n):
    nct = len(city)
    
    newcity=[]
    # Segment in the range n[0]...n[1]
    for j in range( (n[1]-n[0])%nct + 1):
        newcity.append(city[ (j+n[0])%nct ])
    # is followed by segment n[5]...n[2]
    for j in range( (n[2]-n[5])%nct + 1):
        newcity.append(city[ (j+n[5])%nct ])
    # is followed by segment n[3]...n[4]
    for j in range( (n[4]-n[3])%nct + 1):
        newcity.append(city[ (j+n[3])%nct ])
    return newcity

def Plot(city, R, dist):
    # Plot
    Pt = [R[city[i]] for i in range(len(city))]
    Pt += [R[city[0]]]
    Pt = array(Pt)
    title('Total distance='+str(dist))
    plot(Pt[:,0], Pt[:,1], '-o')
    show()

'''
Traveling salesman problem solved using Simulated Annealing.
'''
def pathUpdate(R):
    ncity = 100        # Number of cities to visit
    maxTsteps = 100    # Temperature is lowered not more than maxTsteps
    Tstart = 0.2       # Starting temperature - has to be high enough
    fCool = 0.9        # Factor to multiply temperature at each cooling step
    maxSteps = 100*ncity     # Number of steps at constant temperature
    maxAccepted = 10*ncity   # Number of accepted steps at constant temperature

    Preverse = 0.5      # How often to choose reverse/transpose trial move

    R = array(R)

    # The index table -- the order the cities are visited.
    city = range(ncity)
    # Distance of the travel at the beginning
    dist = TotalDistance(city, R)

    # Stores points of a move
    n = zeros(6, dtype=int)
    nct = len(R) # number of cities
    
    T = Tstart # temperature


    
    for t in range(maxTsteps):  # Over temperature

        accepted = 0
        for i in range(maxSteps): # At each temperature, many Monte Carlo steps
            
            while True: # Will find two random cities sufficiently close by
                # Two cities n[0] and n[1] are choosen at random
                n[0] = int((nct)*rand())     # select one city
                n[1] = int((nct-1)*rand())   # select another city, but not the same
                if (n[1] >= n[0]): n[1] += 1   #
                if (n[1] < n[0]): (n[0],n[1]) = (n[1],n[0]) # swap, because it must be: n[0]<n[1]
                nn = (n[0]+nct -n[1]-1) % nct  # number of cities not on the segment n[0]..n[1]
                if nn>=3: break
        
            # We want to have one index before and one after the two cities
            # The order hence is [n2,n0,n1,n3]
            n[2] = (n[0]-1) % nct  # index before n0  -- see figure in the lecture notes
            n[3] = (n[1]+1) % nct  # index after n2   -- see figure in the lecture notes
            
            if Preverse > rand(): 
                # Here we reverse a segment
                # What would be the cost to reverse the path between city[n[0]]-city[n[1]]?
                de = Distance(R[city[n[2]]],R[city[n[1]]]) + Distance(R[city[n[3]]],R[city[n[0]]]) - Distance(R[city[n[2]]],R[city[n[0]]]) - Distance(R[city[n[3]]],R[city[n[1]]])
                
                if de<0 or exp(-de/T)>rand(): # Metropolis
                    accepted += 1
                    dist += de
                    reverse(city, n)
            else:
                # Here we transpose a segment
                nc = (n[1]+1+ int(rand()*(nn-1)))%nct  # Another point outside n[0],n[1] segment. See picture in lecture nodes!
                n[4] = nc
                n[5] = (nc+1) % nct
        
                # Cost to transpose a segment
                de = -Distance(R[city[n[1]]],R[city[n[3]]]) - Distance(R[city[n[0]]],R[city[n[2]]]) - Distance(R[city[n[4]]],R[city[n[5]]])
                de += Distance(R[city[n[0]]],R[city[n[4]]]) + Distance(R[city[n[1]]],R[city[n[5]]]) + Distance(R[city[n[2]]],R[city[n[3]]])
                
                if de<0 or exp(-de/T)>rand(): # Metropolis
                    accepted += 1
                    dist += de
                    city = transpt(city, n)
                    
            if accepted > maxAccepted: break

            
        print "T=%10.5f , distance= %10.5f , accepted steps= %d" %(T, dist, accepted)
        T *= fCool             # The system is cooled down
        if accepted == 0: 
            break  # If the path does not want to change any more, we can stop

    print city
    print "done"
    return city

# parameters: path of data file, seperator between gene features
def readAll(fname, seperator):
    dataset= []
    with open(fname, 'r') as f:
        for line in f:
            array= [ float(element.strip()) for element in line.split(seperator)]
            dataset.append(array)
    return dataset


def calDistance(path, gene):
    s= 0
    for i, pos in enumerate(path):
        if i==0:
            continue
        s = s+ fabs(gene[i]-gene[i-1])            
    return s

def outputCluster(cluster):
    f= open('clusters.txt', 'w')
    for i in range(len(cluster)):
        if len(cluster[i])>0:
            f.write('%d\n'%(i))
            for j in range(len(cluster[0])):
                
def isStable(prepath, path):
    for i in range(len(path)):
        if prepath[i]!= path[i]:
            return False
    return True                 


def dataGen(nGene, nFeature):
    dataset= [[random() for j in range(nFeature)] for i in range(nGene)]
    return dataset
'''

'''            
def cluster(dataset, nCluster, maxIter):
    size= len(dataset)          #number of genes
    nfeatures= len(dataset[0])  #number of features
    if nCluster> size:
        print 'datasize less than number of clusters'
        exit(1)
        
    cluster= [[] for i in range(nCluster)]
    
    # pick up k genes randomly to init each cluster
    for i in range(nCluster):
        ind= int (random()* len(dataset))
        gene= dataset[ind]
        del dataset[ind]
        cluster[i]= [[gene[j]] for j in len(gene) ]
    
    prepath=[]
    path=[]
    for i in range(nCluster):
        path.append(pathUpdate(cluster[i]))
        prepath=[j for j in range(nfeatures)]  

    iloop= 0              
    while True:        
        minDist= float('inf')
        minCluster= 0  
        for gene in dataset:
            for i in range(nCluster):
                dist= calDistance(path[i], gene)
                if dist< minDist:
                    minDist= dist
                    minCluster= i
            [cluster[minCluster][j].append( gene[j]) for j in len(gene)]       #update cluster
        iloop= iloop+1
        print 'iteration %d\n'%(iloop)
        
        if isStable(path, prepath) or iloop> maxIter:
            break
        else:
            prepath= path
            path[minCluster]= pathUpdate(cluster[minCluster])                
    return cluster
    
'''
dataPath and maxIter are optional
if no dataPath provided, generate random data
'''  
if __name__=='__main__':
    if len(sys.argv)< 2:
        print 'Usage: salesmanKmeans nCluster [dataPath] [maxIter]'
        exit(1)
        
    if len(sys.argv)<3:
        dataset= dataGen(nGene, nFreature)
    else:    
        dataset= readAll(argv[1])       #assume each line represents a gene
        
    if len(sys.argv) ==4:
        maxIter= sys.argv[3]
    else:
        maxIter= 2000       
    
    cluster(dataset)
    

    
    
    
