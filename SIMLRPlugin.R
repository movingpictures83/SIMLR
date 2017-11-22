dyn.load("projsplx_R.so");


# perform the SIMLR clustering algorithm
#"SIMLR" <- function( X, c , no.dim = NA, k = 10, impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {
#"SIMLR" <- function( D_Kernels, c , no.dim = NA, k = 10, impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {
"SIMLR" <- function( D_Kernels, c , pc,  no.dim = NA, k = 10, impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {
    #print("X MATRIX:");
    #print(D_Kernels);
    # set any required parameter to the defaults
    if(is.na(no.dim)) {
        no.dim = c
    }
    
    # check the if.impute parameter
    if(impute == TRUE) {
        X = t(X)
        X_zeros = which(X==0,arr.ind=TRUE)
        if(length(X_zeros)>0) {
            R_zeros = as.vector(X_zeros[,"row"])
            C_zeros = as.vector(X_zeros[,"col"])
            ind = (C_zeros - 1) * nrow(X) + R_zeros
            X[ind] = as.vector(colMeans(X))[C_zeros]
        }
        X = t(X)
    }
    
    # check the normalize parameter
    if(normalize == TRUE) {
        X = t(X)
        X = X - min(as.vector(X))
        X = X / max(as.vector(X))
        C_mean = as.vector(colMeans(X))
        X = apply(X,MARGIN=1,FUN=function(x) return(x-C_mean))
    }
    
    # start the clock to measure the execution time
    ptm = proc.time()
    
    # set some parameters
    NITER = 30
    #num = ncol(X)
    num = ncol(D_Kernels[[1]])
    r = -1
    beta = 0.8
 
    # TMC Commenting Out... Our Plugin Assumes Kernels in Advance   
    #cat("Computing the multiple Kernels.\n")
    # compute the kernels
    #D_Kernels = multiple.kernel(t(X),cores.ratio)
    # set up some parameters
    alphaK = 1 / rep(length(D_Kernels),length(D_Kernels))
    distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
    for (i in 1:length(D_Kernels)) {
        distX = distX + D_Kernels[[i]]
    }
    distX = distX / length(D_Kernels)
    distX = as.matrix(distX)

    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    #print("A")
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    #print("B")
    idx = array(0,c(nrow(distX),ncol(distX)))
    #print("C")
    for(i in 1:nrow(distX)) {
        distX1[i,] = res[[i]]$x
        idx[i,] = res[[i]]$ix
    }
    #print("HELLO")
    A = array(0,c(num,num))
    di = distX1[,2:(k+2)]
    rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
    id = idx[,2:(k+2)]
    
    numerator = (apply(array(0,c(length(di[,k+1]),dim(di)[2])),MARGIN=2,FUN=function(x) {x=di[,k+1]}) - di)
    temp = (k*di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum) + .Machine$double.eps)
    denominator = apply(array(0,c(length(temp),dim(di)[2])),MARGIN=2,FUN=function(x) {x=temp})
    temp = numerator / denominator
    a = apply(array(0,c(length(t(1:num)),dim(di)[2])),MARGIN=2,FUN=function(x) {x=1:num})
    A[cbind(as.vector(a),as.vector(id))] = as.vector(temp)
    if(r<=0) {
        r = mean(rr)
    }
    lambda = max(mean(rr),0)
    A[is.nan(A)] = 0
    A0 = (A + t(A)) / 2
    S0 = max(max(distX)) - distX
    
    cat("Performing network diffiusion.\n")
    
    # perform network diffiusion
    S0 = network.diffusion(S0,k)
    
    # compute dn
    S0 = dn(S0,'ave')
    S = S0
    D0 = diag(apply(S,MARGIN=2,FUN=sum))
    L0 = D0 - S
    
    eig1_res = eig1(L0,c,0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    evs_eig1 = eig1_res$eigval_full
    
    # perform the iterative procedure NITER times
    converge = vector()
    for(iter in 1:NITER) {
        
        cat("Iteration: ",iter,"\n")
        
        distf = L2_distance_1(t(F_eig1),t(F_eig1))
        A = array(0,c(num,num))
        b = idx[,2:dim(idx)[2]]
        a = apply(array(0,c(num,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:num })
        inda = cbind(as.vector(a),as.vector(b))
        ad = (distX[inda]+lambda*distf[inda])/2/r
        dim(ad) = c(num,ncol(b))
        
        # call the c function for the optimization
        c_input = -t(ad)
        c_output = t(ad)
        ad = t(.Call("projsplx_R",c_input,c_output))
        
        A[inda] = as.vector(ad)
        A[is.nan(A)] = 0
        A = (A + t(A)) / 2
        S = (1 - beta) * S + beta * A
        S = network.diffusion(S,k)
        D = diag(apply(S,MARGIN=2,FUN=sum))
        L = D - S
        #print("BEFORE FOLD")
        F_old = F_eig1
        eig1_res = eig1(L,c,0)
        F_eig1 = eig1_res$eigvec
        temp_eig1 = eig1_res$eigval
        ev_eig1 = eig1_res$eigval_full
        evs_eig1 = cbind(evs_eig1,ev_eig1)
        DD = vector()
        #print("AFTER VECTOR")
        #for (i in 1:length(D_Kernels)) {
        for (i in 1:length(D_Kernels)) {
            temp = (.Machine$double.eps+D_Kernels[[i]]) * (S+.Machine$double.eps)
            #print("DD")
            temp = as.matrix(temp)
            DD[i] = mean(apply(temp,MARGIN=2,FUN=sum))
            #print("EE")
        }
        #print("BEFORE ALPHA")
        alphaK0 = umkl(DD)
        alphaK0 = alphaK0 / sum(alphaK0)
        alphaK = (1-beta) * alphaK + beta * alphaK0
        alphaK = alphaK / sum(alphaK)
        fn1 = sum(ev_eig1[1:c])
        fn2 = sum(ev_eig1[1:(c+1)])
        converge[iter] = fn2 - fn1
        if (iter<10) {
            if (ev_eig1[length(ev_eig1)] > 0.000001) {
                lambda = 1.5 * lambda
                r = r / 1.01
            }
        }
        else {
            if(converge[iter]>converge[iter-1]) {
                S = S_old
                if(converge[iter-1] > 0.2) {
                    warning('Maybe you should set a larger value of c.')
                }
                break
            }
        }
        S_old = S
        
        # compute Kbeta
        distX = D_Kernels[[1]] * alphaK[1]
        #print("RUNNING KERNELS")
        for (i in 2:length(D_Kernels)) {
            distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
        }
        #print("DONE RUNNING KERNELS")
        distX = as.matrix(distX)
        # sort distX for rows
        res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
        #print("AFTER RES")
        distX1 = array(0,c(nrow(distX),ncol(distX)))
        idx = array(0,c(nrow(distX),ncol(distX)))
        for(i in 1:nrow(distX)) {
            distX1[i,] = res[[i]]$x
            idx[i,] = res[[i]]$ix
        }
        
    }
    LF = F_eig1
    D = diag(apply(S,MARGIN=2,FUN=sum))
    L = D - S
    
    # compute the eigenvalues and eigenvectors of P
    eigen_L = eigen(L)
    U = eigen_L$vectors
    D = eigen_L$values
    
    if (length(no.dim)==1) {
        U_index = seq(ncol(U),(ncol(U)-no.dim+1))
        F_last = tsne(S,k=no.dim,initial_config=U[,U_index])
    }
    else {
        F_last = list()
        for (i in 1:length(no.dim)) {
            U_index = seq(ncol(U),(ncol(U)-no.dim[i]+1))
            F_last[i] = tsne(S,k=no.dim[i],initial_config=U[,U_index])
        }
    }
    
    # compute the execution time
    execution.time = proc.time() - ptm
    
    cat("Performing Kmeans.\n")
    y = kmeans(F_last,c,nstart=200)
    #y = kmeans(pc, c, nstart=200)    

    ydata = tsne(S)
    
    # create the structure with the results
    results = list()
    results[["y"]] = y
    results[["S"]] = S
    results[["F"]] = F_last
    results[["ydata"]] = ydata
    results[["alphaK"]] = alphaK
    results[["execution.time"]] = execution.time
    results[["converge"]] = converge
    results[["LF"]] = LF
    #print("Y: ");
    #print(y);
    return(results)
    
}


# perform network diffusion of K steps over the network A
"network.diffusion" <- function( A, K ) {

    # set the values of the diagonal of A to 0
    diag(A) = 0

    # compute the sign matrix of A
    sign_A = A
    sign_A[which(A>0,arr.ind=TRUE)] = 1
    sign_A[which(A<0,arr.ind=TRUE)] = -1

    # compute the dominate set for A and K
    P = dominate.set(abs(A),min(K,nrow(A)-1)) * sign_A

    # sum the absolute value of each row of P
    DD = apply(abs(P),MARGIN=1,FUN=sum)

    # set DD+1 to the diagonal of P
    diag(P) = DD + 1

    # compute the transition field of P
    P = transition.fields(P)

    # compute the eigenvalues and eigenvectors of P
    eigen_P = eigen(P)
    U = eigen_P$vectors
    D = eigen_P$values

    # set to d the real part of the diagonal of D
    d = Re(D + .Machine$double.eps)

    # perform the diffusion
    alpha = 0.8
    beta = 2
    d = ((1-alpha)*d)/(1-alpha*d^beta)

    # set to D the real part of the diagonal of d
    D = array(0,c(length(Re(d)),length(Re(d))))
    diag(D) = Re(d)

    # finally compute W
    W = U %*% D %*% t(U)
    diagonal_matrix = array(0,c(nrow(W),ncol(W)))
    diag(diagonal_matrix) = 1
    W = (W * (1-diagonal_matrix)) / apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=(1-diag(W))})
    diag(D) = diag(D)[length(diag(D)):1]
    W = diag(DD) %*% W
    W = (W + t(W)) / 2

    W[which(W<0,arr.ind=TRUE)] = 0

    return(W)

}

# compute the dominate set for the matrix aff.matrix and NR.OF.KNN
"dominate.set" <- function( aff.matrix, NR.OF.KNN ) {

    # create the structure to save the results
    PNN.matrix = array(0,c(nrow(aff.matrix),ncol(aff.matrix)))

    # sort each row of aff.matrix in descending order and saves the sorted 
    # array and a collection of vectors with the original indices
    res.sort = apply(t(aff.matrix),MARGIN=2,FUN=function(x) {return(sort(x, decreasing = TRUE, index.return = TRUE))})
    sorted.aff.matrix = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$x) }))
    sorted.indices = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$ix) }))

    # get the first NR.OF.KNN columns of the sorted array
    res = sorted.aff.matrix[,1:NR.OF.KNN]

    # create a matrix of NR.OF.KNN columns by binding vectors of 
    # integers from 1 to the number of rows/columns of aff.matrix
    inds = array(0,c(nrow(aff.matrix),NR.OF.KNN))
    inds = apply(inds,MARGIN=2,FUN=function(x) {x=1:nrow(aff.matrix)})

    # get the first NR.OF.KNN columns of the indices of aff.matrix
    loc = sorted.indices[,1:NR.OF.KNN]

    # assign to PNN.matrix the sorted indices
    PNN.matrix[(as.vector(loc)-1)*nrow(aff.matrix)+as.vector(inds)] = as.vector(res)

    # compute the final results and return them
    PNN.matrix = (PNN.matrix + t(PNN.matrix))/2

    return(PNN.matrix)

}


# compute the transition field of the given matrix
"transition.fields" <- function( W ) {

    # get any index of columns with all 0s
    zero.index = which(apply(W,MARGIN=1,FUN=sum)==0)

    # compute the transition fields
    W = dn(W,'ave')

    w = sqrt(apply(abs(W),MARGIN=2,FUN=sum)+.Machine$double.eps)
    W = W / t(apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=w}))
    W = W %*% t(W)

    # set to 0 the elements of zero.index
    W[zero.index,] = 0
    W[,zero.index] = 0

    return(W)

}


# normalizes a symmetric kernel
"dn" = function( w, type ) {

    # compute the sum of any column
    D = apply(w,MARGIN=2,FUN=sum)

    # type "ave" returns D^-1*W
    if(type=="ave") {
        D = 1 / D
        D_temp = matrix(0,nrow=length(D),ncol=length(D))
        D_temp[cbind(1:length(D),1:length(D))] = D
        D = D_temp
        wn = D %*% w
    }
    # type "gph" returns D^-1/2*W*D^-1/2
    else if(type=="gph") {
        D = 1 / sqrt(D)
        D_temp = matrix(0,nrow=length(D),ncol=length(D))
        D_temp[cbind(1:length(D),1:length(D))] = D
        D = D_temp
        wn = D %*% (w %*% D)
    }
    else {
        stop("Invalid type!")
    }

    return(wn)

}

# normalizes a symmetric kernel
"dn" = function( w, type ) {

    # compute the sum of any column
    D = apply(w,MARGIN=2,FUN=sum)

    # type "ave" returns D^-1*W
    if(type=="ave") {
        D = 1 / D
        D_temp = matrix(0,nrow=length(D),ncol=length(D))
        D_temp[cbind(1:length(D),1:length(D))] = D
        D = D_temp
        wn = D %*% w
    }
    # type "gph" returns D^-1/2*W*D^-1/2
    else if(type=="gph") {
        D = 1 / sqrt(D)
        D_temp = matrix(0,nrow=length(D),ncol=length(D))
        D_temp[cbind(1:length(D),1:length(D))] = D
        D = D_temp
        wn = D %*% (w %*% D)
    }
    else {
        stop("Invalid type!")
    }

    return(wn)

}


# compute the eigenvalues and eigenvectors
"eig1" <- function( A, c = NA, isMax = NA, isSym = NA ) {

    # set the needed parameters
    if(is.na(c)) {
        c = dim(A)[1]
    }
    if(c>dim(A)[1]) {
        c = dim(A)[1]
    }
    if(is.na(isMax)) {
        isMax = 1
    }
    if(is.na(isSym)) {
        isSym = 1
    }

    # compute the eigenvalues and eigenvectors of A
    if(isSym==1) {
        eigen_A = eigen(A,symmetric=TRUE)
    }
    else {
        eigen_A = eigen(A)
    }
    v = eigen_A$vectors
    d = eigen_A$values

    # sort the eigenvectors
    if(isMax == 0) {
        eigen_A_sorted = sort(d,index.return=TRUE)
    }
    else {
        eigen_A_sorted = sort(d,decreasing=TRUE,index.return=TRUE)
    }
    d1 = eigen_A_sorted$x
    idx = eigen_A_sorted$ix
    idx1 = idx[1:c]

    # compute the results
    eigval = d[idx1]
    eigvec = Re(v[,idx1])
    eigval_full = d[idx]

    return(list(eigval=eigval,eigvec=eigvec,eigval_full=eigval_full))

}


# compute the L2 distance
"L2_distance_1" <- function( a, b ) {

    if(dim(a)[1] == 1) {
        a = rbind(a,rep(0,dim(a)[2]))
        b = rbind(b,rep(0,dim(b)[2]))
    }

    aa = apply(a*a,MARGIN=2,FUN=sum)
    bb = apply(b*b,MARGIN=2,FUN=sum)
    ab = t(a) %*% b
    d1 = apply(array(0,c(length(t(aa)),length(bb))),MARGIN=2,FUN=function(x){ x = t(aa) })
    d2 = t(apply(array(0,c(length(t(bb)),length(aa))),MARGIN=2,FUN=function(x){ x = t(bb) }))
    d = d1 + d2 - 2 * ab
    d = Re(d)
    d = matrix(mapply(d,FUN=function(x) { return(max(max(x),0)) }),nrow=nrow(d),ncol=ncol(d))
    d_eye = array(1,dim(d))
    diag(d_eye) = 0
    d = d * d_eye

    return(d)

}


# umkl function
"umkl" = function( D, beta = NA ) {

    # set some parameters
    if(is.na(beta)) {
        beta = 1 / length(D)
    }
    tol = 1e-4
    u = 20
    logU = log(u)

    # compute Hbeta
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P

    betamin = -Inf
    betamax = Inf
    # evaluate whether the perplexity is within tolerance
    Hdiff = H - logU
    tries = 0
    while (abs(Hdiff) > tol && tries < 30) {
        #if not, increase or decrease precision
        if (Hdiff > 0) {
            betamin = beta
            if(abs(betamax)==Inf) {
                beta = beta * 2
            }
            else {
                beta = (beta + betamax) / 2
            }
        }
        else {
            betamax = beta
            if(abs(betamin)==Inf) {
                beta = beta * 2
            }
            else {
                beta = (beta + betamin) / 2
            }
        }
        # compute the new values
        res_hbeta = Hbeta(D, beta)
        H = res_hbeta$H
        thisP = res_hbeta$P
        Hdiff = H - logU
        tries = tries + 1
    }

    return(thisP)

}

"Hbeta" = function( D, beta ) {

    D = (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
    P = exp(-D * beta)
    sumP = sum(P)
    H = log(sumP) + beta * sum(D * P) / sumP
    P = P / sumP

    return(list(H=H,P=P))

}


"tsne" <- function(X, initial_config = NULL, k = 2, max_iter = 1000, min_cost = 0, epoch = 100) {

    cat("Performing t-SNE.\n")

    momentum = 0.8
    final_momentum = 0.8
    mom_switch_iter = 250

    epsilon = 500
    min_gain = 0.01
    initial_P_gain = 4

    n = nrow(X)

    eps = .Machine$double.eps

    if (!is.null(initial_config) && is.matrix(initial_config)) {
        if (nrow(initial_config) != n | ncol(initial_config) != k) {
            stop('initial_config argument does not match necessary configuration for X')
        }
        ydata = initial_config
        initial_P_gain = 1

    }
    else {
        ydata = matrix(rnorm(k * n),n)
    }

    P = X
    P = 0.5 * (P + t(P))

    P[P < eps]<-eps
    P = P/sum(P)

    P = P * initial_P_gain
    grads = matrix(0,nrow(ydata),ncol(ydata))
    incs = matrix(0,nrow(ydata),ncol(ydata))
    gains = matrix(1,nrow(ydata),ncol(ydata))

    for (iter in 1:max_iter) {

        if (iter %% epoch == 0) {
            cost = sum(apply(P * log((P+eps)/(Q+eps)),1,sum))

            cat("Epoch: Iteration #",iter," error is: ",cost,"\n")

            if (cost < min_cost) {
                break
            }
        }

        sum_ydata = apply((ydata^2),1,sum)
        num =  1/(1 + sum_ydata + sweep(-2*ydata %*% t(ydata),2,-t(sum_ydata)))
        diag(num) = 0
        Q = num / sum(num)

        if (any(is.nan(num))) {
            message ('NaN in grad. descent')
        }

        Q[Q < eps] = eps

        stiffnesses = (P-Q) * num
        grads = 4 * (diag(apply(stiffnesses,2,sum)) - stiffnesses) %*% ydata
        gains = (gains + .2) * abs(sign(grads) != sign(incs)) + gains * .8 * abs(sign(grads) == sign(incs))
        gains[gains < min_gain] = min_gain

        incs = momentum * incs - epsilon * (gains * grads)
        ydata = ydata + incs
        ydata = sweep(ydata,2,apply(ydata,2,mean))

        # we are constraining the ydata
        ydata[ydata < -100] = -100
        ydata[ydata > 100] = 100

        if (iter == mom_switch_iter) {
            momentum = final_momentum
        }

        if (iter == 100 && is.null(initial_config)) {
            P = P/4
        }

    }

    return(ydata)

}






p_value <- 0.01;
libs <- c("Hmisc");
lapply(libs, require, character.only=T);


input <- function(inputfile) {
  #X <<- array(dim=4);
  #X[1] <<- as.matrix(read.csv(paste(inputfile, ".pearson.csv", sep=""), header = TRUE));
  #X[2] <<- as.matrix(read.csv(paste(inputfile, ".spearman.csv", sep=""), header = TRUE));
  #X[3] <<- as.matrix(read.csv(paste(inputfile, ".dcov.csv", sep=""), header = TRUE));
  #X[4] <<- as.matrix(read.csv(paste(inputfile, ".MIC.csv", sep=""), header = TRUE));
  pc <<- read.csv(paste(inputfile, ".pearson.csv", sep=""), header = TRUE);
 
  X <<- list(as.matrix(pc), as.matrix(pc)
             #as.matrix(read.csv(paste(inputfile, ".spearman.csv", sep=""), header = TRUE)),
             #as.matrix(read.csv(paste(inputfile, ".dcov.csv", sep=""), header = TRUE)),
             #as.matrix(read.csv(paste(inputfile, ".MIC.csv", sep=""), header = TRUE))
            );
  for (k in 1:length(X)) {
     for (i in 1:length(pc)) {
        for (j in 1:length(pc)) {
           print("BEFORE");
           print(X[[k]][i, j]);
           X[[k]][i, j] = 1 - X[[k]][i, j]
           print("AFTER");
           print(X[[k]][i, j]);
        }
     }
  }
}


run <- function() {
  cn <<- colnames(pc);
  #cn <<- cn[2:length(cn)];
  #D <- 0.25*D_pearson + 0.25*D_spearman + 0.25*D_dcov + 0.25*D_MIC;
  nclust <<- 23
  results <<- SIMLR(X, nclust, as.matrix(pc));  # Just for now, will vary cluster count later
}

output <- function(outputfile) {
   ##print("xxxxxxxxxxxxxxxxxxxxxxxxxx");
   #print("CN");
   #print(cn);
   #print("Y:");
   #print(results$y);
   #print("CLUSTER:");
   #print(results$y$cluster);
   print("LENGTH");
   print(length(results$y$cluster));
   #print("HELLO");
   #fileConn<-file(outputfile);
   #cat("hello");
   #cat("\n");
   #cat("world");
   #write("hello\n", fileConn);
   #write("world\n", fileConn);
   outNOA <- paste(outputfile,".noa", sep="");
   print("NCLUST");
   print(nclust);
   print(outNOA);
   myClust <- vector("list", nclust);
   print("Z");
   print("RESULTSYCLUSTER");
   print(results$y$cluster);
   sink(outNOA);
   #print("Y");
   cat("Name");
   #print("X");
   cat("\t");
   #print("W");
   cat("Cluster");
   #print("V");
   cat("\n");
   #print("HELLO");
   #print(length(results$y$cluster));
   for (i in 1:length(results$y$cluster)) {
      #print("A");
      cat(noquote(names(pc))[i]);
      #print("B");
      cat("\t");
      #print("C");
      cat(strtoi(results$y$cluster[i]));
      #print("D");
      cat("\n");
      #print("BEFORE CLUSTER");
      myClust[[strtoi(results$y$cluster[i])]] = c(myClust[[strtoi(results$y$cluster[i])]], names(pc)[i]);
      #print("INTERNAL CLUSTER");
      #print(strtoi(results$y$cluster[i]));
      #print(myClust[[strtoi(results$y$cluster[i])]]);
   }
   sink();
   #print("DOING CSV....");
   outCSV <- paste(outputfile,".csv",sep="");
   sink(outCSV);
   for (i in 1:nclust) {
      cat("\"\",");
      cat("\t");
      cat("\"x\"");
      cat("\n");
      for (j in 1:length(myClust[[i]])) {
         cat("\"");
         cat(j);
         cat("\"");
         cat(",");
         cat(myClust[[i]][j]);
         cat("\n");
      }
   }
   sink();
   #close(fileConn);
   #print("FIRST ROW");
   #print(names(pc)[1]);
   #print("xxxxxxxxxxxxxxxxxxxxxxxxxx");
   #write.table(results$S, file=outputfile, sep=",", append=FALSE, row.names=unlist(cn), col.names=unlist(cn), na="");
}


