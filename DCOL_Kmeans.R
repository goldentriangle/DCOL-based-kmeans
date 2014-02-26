#implementation of DCOL-based kmeans clustering
#contact: kai.wang.magic@gmail.com  Kai

library(TSP)
library(doSNOW)
#library(clues)
library(dynamicTreeCut)
library(fossil)
library(phyclust, quiet = TRUE)
# to use function adjustedRand() for quality assessment

gene.specific.null<-function(array, B=500)
{
	null.mat<-matrix(0, nrow=nrow(array), ncol=B)
	l<-ncol(array)
	d.array<-array[,1:(l-1)]
	for(i in 1:B)
	{
		this.order<-sample(l, l, replace=FALSE)
		for(j in 1:(l-1)) d.array[,j]<-abs(array[,this.order[j+1]]-array[,this.order[j]])
		null.mat[,i]<-apply(d.array, 1, sum)
	}
	r<-cbind(apply(null.mat, 1, mean), apply(null.mat, 1, sd))
	return(r)
}

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
        #this.size<-aver.grp.size
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
	r$data<-a
	r$grps<-g2
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
          

data.cluster<-function( dataset, nCluster, grps, maxIter =10){  
    #dataset= x#$data
    #grps = x$grps
    
    null.mat<- gene.specific.null(dataset)  #mean and sd for each gene
    
    dims= dim(dataset)  
    cat('dim:',dims,'\n')        
    ngene= dims[1] 
    nsample= dims[2]
    
    if (nCluster > ngene){
        sprintf('datasize %d smaller than the amount of clusters %d\n',nCluster, size)
        q()
    }
    
    dataCluster<- list()
    indCluster<- list()
    # pick up one gene for each cluster for initialization randomly
    ind= sample(1: ngene, nCluster)
    
    prepath= list()
    path= list()
    
    cnt=1
    for(i in ind){
        #print(matrix(dataset[i,], nrow= nsample))
        d= dist( matrix(dataset[i,], nrow= nsample))
        #print (d)
        tsp= TSP(d)
        tour= solve_TSP(tsp)   
        #print(as.integer(tour))        
        path[[cnt]]= as.integer(tour)
        cnt= cnt+1
    }    
    #print(prepath)
    #print(path)
    
    p_max= 0.2            # initial p-value threshold
    p_min= 0.05           # more strict p-value
    p_val= p_max
    
    iloop= 0 
    unclustered=c()             
    while (TRUE){  
        dataCluster= list()
        indCluster= list()          
        isEmpty= rep(TRUE, nCluster)
        unclustered=c()
        
        iloop= iloop+1
        if (iloop%%5 == 1){
            cat('------clustering iteration ',iloop,'-------\n')
        }
        
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
            
            #test p-value
            if(pnorm(minDist, mean=null.mat[pnt_i,1], sd=null.mat[pnt_i,2], lower.tail=TRUE)>p_val){
                unclustered= c(unclustered, pnt_i)           
                next                           
            }
            
            #print(isEmpty)
            if( isEmpty[minCluster]){
                #print('branch 1')
                dataCluster[[minCluster]]= rbind(NULL, dataset[pnt_i,])
                indCluster[[minCluster]]= pnt_i
                isEmpty[minCluster]= FALSE
            }else{
                #print('branch 2')
                dataCluster[[minCluster]]= rbind(dataCluster[[minCluster]], dataset[pnt_i,])
                #print(indCluster[[minCluster]])
                indCluster[[minCluster]]= c(indCluster[[minCluster]],pnt_i)
                #print(indCluster[[minCluster]])
            }     
        }
        
        p_val <- p_val- (p_max-p_min)/maxIter
      
        prepath= path
        for(i in 1:nCluster){
            if(isEmpty[i]== FALSE){
                d= dist(t(dataCluster[[i]]))
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
    cat('-----------------simulated results----------------\n')
    r<-new("list")
    r$cluster= matrix(0, ngene, 1)
    
    for(i in 1:nCluster){
        cat('cluster ',i, '\n')
        if(isEmpty[i]== FALSE){
            a<- grps[indCluster[[i]] ]
            #print(a)
            r$cluster[indCluster[[i]] ] = i
            
            cat('the most frequent class falling into this group is', names(which.max(table(a))),'\n')
            cat('majority/total is',max(table(a)),'/',length(indCluster[[i]]),'=',max(table(a))/length(indCluster[[i]]),'\n\n\n')
        }
        
    }
    cat('\nUnclustered:', length(unclustered),':\n')
    #print(grps[unclustered])
    print('ret')
    return (r)
}

# find the best K of k-means clustering by Gap Statistic
findK <- function(dataset){
    source('clusGap.R')
    n.genes<-1000
    n.samples<-c(100)
    n.grps<-5
    aver.grp.size<-200
    n.fun.types<-c(4)
    epsilon=c(0.2)
    n.depend= 2

    dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)

    if (0 %in% unique(dataset$grps)){
        ngrp= length(unique(dataset$grps))-1
    }else{
        ngrp= length(unique(dataset$grps))
    }
    print(ngrp)

    gskmn <- clusGap(dataset$data, FUN = data.cluster, grps=dataset$grps, K.max = ngrp+3, B = 60)
    gskmn
    plot(gskmn, main = "clusGap(., FUN = kmeans, n.start=20, B= 60)")
}

test <- function(){
	sink("log.txt", split="TRUE")
	f= file(sprintf('result%s.txt', Sys.time()),'a')
	write.table(sprintf("\n\n\nstart running at %s .....\n",Sys.time()), f)
	for (n.depend in 1:1){
		for (i in seq(2,7)){             #epsilon ranges from 0.2 to 0.7
		n.samples<-c(100)
	    n.grps<-10                  #total number of genes equal to n.grps * aver.grp.size
	    aver.grp.size<-100
	    n.genes<- n.grps * aver.grp.size
	    #n.genes<- (1+n.grps) * aver.grp.size
	    n.fun.types<-c(4)
	    epsilon=c(0.1*i)
	     
	
	    dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
		cat('epsilon: ',epsilon, '\n')
		cat('depend: ',n.depend, '\n')
		writeLines(sprintf("epsilon %f, depend %f\n", epsilon, n.depend),f)
		flush(f)
		
		dr<- c()
		kr<- c()
		hr<- c()
		#change number of iterations for each parameter setting here
		for(j in 1:2){
			write(sprintf("testing iteration %d", j), stderr())

			if (0 %in% unique(dataset$grps)){
		        ngrp= length(unique(dataset$grps))-1
		    }else{
		        ngrp= length(unique(dataset$grps))
		    }
		    print(ngrp)
		    #cat ("truth\n",dataset$grps,'\n')
            
            if(TRUE){
			    results= data.cluster(dataset$data, ngrp, dataset$grps, 10 )$cluster
			    #cat ('DCOL\n', '\n') #results,
			    write("DCOL", stderr())
			    print(table (dataset$grps, results))
			    print(RRand(dataset$grps +1 , results+1))
			    dr<- append(dr, RRand(dataset$grps +1 , results+1)$adjRand)
	            print(dr[length(dr)])
				
							
	            #cat("kmeans\n")
			    write("kmeans", stderr())
			    kresult= kmeans(dataset$data, ngrp)
			    print(table (dataset$grps, kresult$cluster))
			    print(RRand(dataset$grps +1 , kresult$cluster))
	        
			    kr<- append(kr, RRand(dataset$grps +1 , kresult$cluster)$adjRand)
			    print(kr[length(kr)])
			
			}else{
                source('GDHC.R')
                #cat("GDHC\n")
                write("GDHC", stderr())
                hc= gdhc(dataset$data)
                hgrps= cutreeDynamic(hc, minClusterSize=50)
                #hgrps<-cutree(hc, h=log10(0.01))
                
                print(hgrps)
                print(table (dataset$grps, hgrps))
                print(RRand(dataset$grps + 1, hgrps+1))
                hr<- append(hr, RRand(dataset$grps + 1, hgrps+1)$adjRand)
                print(hr[length(hr)])
            }
		}
        if(TRUE){
            #last column is the mean of all runs
            write.table(cbind( rbind(dr, kr), rowMeans(rbind(dr, kr))),f)
            print(cbind( rbind(dr, kr), rowMeans(rbind(dr, kr))))
        }else{
            # write.table(cbind( rbind(dr, kr), rowMeans(rbind(dr, kr))),f)
            # print(cbind( rbind(dr, kr), rowMeans(rbind(dr, kr))))
            write.table(cbind( rbind(hr), rowMeans(rbind(hr))),f)
            print(cbind( rbind(hr), rowMeans(rbind(hr))))
        }
        
	   	flush(f)
	}
	}
	
	
	close(f)
    sink()

}


test()
proc.time()
