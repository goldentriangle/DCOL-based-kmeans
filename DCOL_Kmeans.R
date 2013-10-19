#implementation of DCOL-based kmeans clustering
#contact: kai.wang.magic@gmail.com  Kai

library(TSP)
library(doSNOW)
library(clues)
# to use function adjustedRand() for quality assessment


data.gen<-function(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
{
	link<-function(x, type)
	{
		x<-(x-mean(x))/sd(x)
		if(type == 1) return(x)
		if(type == 2) return(sin(2*x))
		if(type == 3) return(x^2)
		if(type == 4) return(abs(x))
		if(type == 5) return(x^3)
		if(type == 6) return(atan(4*x))
	}
	
	a<-matrix(rnorm(n.genes*n.samples),ncol=n.samples)
	curr.count<-0
	g<-new("list")
	for(i in 1:n.grps)
	{
		this.size<-rpois(1, aver.grp.size)
		if(this.size < 2) this.size<-2
		
		this.mat<-matrix(0, nrow=this.size, ncol=n.samples)
		this.mat[1,]<-rnorm(n.samples)
		for(j in 2:this.size)
		{
			if(n.depend==0)
			{
				this.basis<-c(1, rep(0,j-2))
			}else{
#				this.basis<-sample(c(1,0), j-1, replace=T, prob=c(min(1, n.depend/(j-1)), 1-min(1, n.depend/(j-1))))
				if(j-1 <= n.depend) 
				{
					this.basis<-rep(1, j-1)
				}else{
					this.basis<-sample(c(rep(1, n.depend), rep(0,j-1-n.depend)), j-1, replace=F)
				}
				
			}
			if(sum(this.basis) > 0)
			{
				x<-rep(0,n.samples)
				for(k in which(this.basis == 1))
				{
					x<-x+link(this.mat[k,], sample(n.fun.types,1))*runif(1,min=-1,max=1)
				}
#				x[x>quantile(x, 0.95)]<-quantile(x, 0.95)
#				x[x<quantile(x, 0.05)]<-quantile(x, 0.05)
				this.mat[j,]<-x
				this.mat[j,]<-(this.mat[j,]-mean(this.mat[j,]))/sd(this.mat[j,])
			}else{
				this.mat[j,]<-rnorm(n.samples)
			}
		}
		if(n.depend == 0)
		{
			this.mat[1,]<-link(this.mat[1,], sample(n.fun.types,1))
			this.mat[1,]<-(this.mat[1,]-mean(this.mat[1,]))/sd(this.mat[1,])
		}
		
		if(curr.count+this.size <= n.genes)
		{
			a[(curr.count+1):(curr.count+this.size),]<-this.mat
			g[[length(g)+1]]<-(curr.count+1):(curr.count+this.size)
		}
		curr.count<-curr.count+this.size		
	}
	a<-a+matrix(rnorm(n.genes*n.samples, sd=epsilon),ncol=n.samples)
	
	g2<-rep(0, nrow(a))
	for(i in 1:length(g)) g2[g[[i]]]<-i
	
	r<-new("list")
	r$data<-a[which(g2>0),]
	r$grps<-g2[which(g2>0)]
	return(r)
}


calDistance<- function(path, gene){
    path= strtoi(path)
    len= length(path)
    #print(len)
    if(len<2){
        print('number of samples < 2')
        q()           
    }
    d= 0
    for(i in 2:len){
        d = d+ abs(gene[path[i]]-gene[path[i-1]])            
    }
    return( d)
}


outputCluster<- function(cluster, path, indCluster){
   
    f= file('clusters.txt')
    for(i in 1:length(cluster)){
        sprintf('cluster %d contains %d members\n', i, length(cluster[[i]]) )
        #cat(path[[i]],'\n')
        #write('%d\n',i)
        
        nrow= dim(cluster[[i]])[1]
        print()          
    }            
}
          

data.cluster<-function(dataset,grps, nCluster, maxIter){    
    dims= dim(dataset)  
    cat('dim:',dims)        
    ngene= dims[1] 
    nsample= dims[2]
    
    if (nCluster > ngene){
        sprintf('datasize %d smaller than the amount of clusters %d\n',nCluster, size)
        q()
    }
    
    cluster<- list()
    indCluster<- list()
    # pick up one gene for each cluster for initialization randomly
    ind= sample(1: ngene, nCluster)
    
    prepath= list()
    path= list()
    
    cnt=1
    for(i in ind){
        print(matrix(dataset[i,], nrow= nsample))
        d= dist( matrix(dataset[i,], nrow= nsample))
        #print (d)
        tsp= TSP(d)
        tour= solve_TSP(tsp)   
        print(as.integer(tour))        
        path[[cnt]]= as.integer(tour)
        cnt= cnt+1
    }    
    print(prepath)
    print(path)
    
    iloop= 1              
    while (TRUE){  
        cluster= list()
        indCluster= list()          
        isEmpty= rep(TRUE, nCluster)
        for(pnt_i in 1:ngene){
            minDist= Inf
            minCluster= 1
           
            #print 'point %d'%(pnt_i)
            for (i in 1:nCluster){
                dist= calDistance(path[[i]], dataset[pnt_i,])
                #cat( dist, minDist, '\n')
                if( dist< minDist){
                    minDist= dist
                    minCluster= i
                }                    
            }     
            #print(isEmpty)
            if( isEmpty[minCluster]){
                #print('branch 1')
                cluster[[minCluster]]= rbind(NULL, dataset[pnt_i,])
                indCluster[[minCluster]]= pnt_i
                isEmpty[minCluster]= FALSE
            }else{
                #print('branch 2')
                cluster[[minCluster]]= rbind(cluster[[minCluster]], dataset[pnt_i,])
                #print(indCluster[[minCluster]])
                indCluster[[minCluster]]= c(indCluster[[minCluster]],pnt_i)
                #print(indCluster[[minCluster]])
            }     
        }
        
        iloop= iloop+1
        if (iloop%%5 == 1){
            cat('------iteration ',iloop,'-------\n')
        }
        prepath= path
        for(i in 1:nCluster){
            if(isEmpty[i]== FALSE){
                d= dist(t(cluster[[i]]))
                tsp= TSP(d)
                tour= solve_TSP(tsp)           
                path[[i]]= as.integer(tour)
            }
        }
          
        if( identical(path, prepath) || iloop> maxIter){
            break
        }           
    }                    
    #return [cluster, path, indCluster]
    #print(path)
    #print(indCluster)
    for(i in 1:nCluster){
        if(isEmpty[i]== FALSE){
            print(grps[indCluster[[i]] ])
        }
        
    }
}

n.genes<-500
n.samples<-c(100)
n.grps<-5
aver.grp.size<-100
n.fun.types<-c(4)
epsilon=c(0.2)
n.depend= 2

dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)

print(length(unique(dataset$grps)))
data.cluster(dataset$data, dataset$grps, length(unique(dataset$grps)),100 )
