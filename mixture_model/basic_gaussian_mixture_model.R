library(pdfCluster)
data("oliveoil")


run_basic_gaussian_mixture <- function(x,g,max_iter = 3, k=NULL,encode=FALSE, seed=1){

  x <- as.matrix(x)
  
  if (encode){
    config <- gcode::extract_config(F)
    config$init <- list(alpha="rnorm",beta="rnorm")
    config$i_dim <- k
    config$j_dim <- k
    config$verbose <- F
    config$tol <- 1e-5
    config$max_iter <- 1000
    config$n.cores <- 8
    config$learn_rate <- 0.0001
    config$batch_size <- k
    config$seed <- seed
    
    join <- gcode::extract_join_framework(F)
    join$complete <- lapply(join$complete,function(X){c(1)})
    join$covariance <- c(1)
    
    gcode.model <- gcode::gcode(data_list = list(t(x)%*%x), config = config, join = join)
    feature_encode <- (gcode.model$main.parameters$beta_sample[[1]])
    
  }
  
  calculate_posterior <- function(x,mu,Sigma,prob){
    posterior <- c()
    
    for (i in c(1:g)){
      posterior[[i]] <- log(prob[i])+mvtnorm::dmvnorm(x,mean = mu[[i]],sigma = Sigma[[i]], log = T)
    }
    
    pys <- do.call('cbind',lapply(c(1:g),function(X){
      main_log_prob <- posterior[[X]] - Reduce("+",posterior)
    }))
    
    pys_max <- apply(pys, 1, max)
    pys <- sweep(pys, 1, pys_max, '-')
    pys <- exp(pys)
    ps.y <- sweep(pys, 1, rowSums(pys), '/')
    
    posterior <- lapply(c(1:3),function(X){ps.y[,X]})
    return(posterior)
  }
  
  
  
  
  
  
  
  
  
  init_cluster <- kmeans(x,g)
  
  prob <- c(table(init_cluster$cluster)/dim(x)[1])
  
  
  mu <- lapply(c(1:g),function(X){
    colMeans(x[init_cluster$cluster==X,])
  })
  
  Sigma <- lapply(c(1:g),function(X){
    cov(x[init_cluster$cluster==X,])
  })
  
  if(encode){
    x_encode <- x%*%feature_encode
    
    mu_encode <- mu
    Sigma_encode <- Sigma
    
    for (i in c(1:g)){
      
      mu_encode[[i]] <- t(feature_encode)%*%array(mu[[i]],dim=c(length(mu[[i]]),1))
      Sigma_encode[[i]] <- t(feature_encode)%*%Sigma[[i]]%*%(feature_encode)
      
    } 
    
  }
  
  for (i in 1:max_iter){
    if (encode){
      posterior <- calculate_posterior(x_encode,mu_encode,Sigma_encode,prob)
      
    } else {
      posterior <- calculate_posterior(x,mu,Sigma,prob)
      
    }
    
    prob <- do.call('c',lapply(posterior,"mean"))
    
    for (i in c(1:g)){
      mu[[i]] <- colMeans((posterior[[i]]*x))/sum(Reduce("+",posterior))
      Sigma[[i]] <- (t(sqrt(posterior[[i]])*x-mu[[i]])%*%(sqrt(posterior[[i]])*x-mu[[i]]))/sum(Reduce("+",posterior))
      
      if(encode){
        mu_encode[[i]] <- t(feature_encode)%*%array(mu[[i]],dim=c(length(mu[[i]]),1))
        Sigma_encode[[i]] <- t(feature_encode)%*%Sigma[[i]]%*%(feature_encode)
        
      }
      
    }
    
  }
  
  
  return(list(mu=mu,Sigma=Sigma,prob=prob,posterior=posterior))
}


p_frame <- c()
for (i in c(1:1000)){
  print(i)
  p <- c()
  basic_mixture <- run_basic_gaussian_mixture(x=oliveoil[,-c(1,2)],g=3,max_iter = 3 ,seed = i)
  p[["basic"]] <- mclust::adjustedRandIndex(apply(do.call('cbind',basic_mixture$posterior),1,function(X){which(X==max(X))}), y = oliveoil[,1])
  
  for (k in c(2:8)){
    
    encoded_mixture <- run_basic_gaussian_mixture(x=oliveoil[,-c(1,2)], g=3, max_iter = 3, k = k,encode = T, seed = i)
    p[[k]] <- mclust::adjustedRandIndex(apply(do.call('cbind',encoded_mixture$posterior),1,function(X){which(X==max(X))}), y = oliveoil[,1])
    
  }
  names(p) <- c("basic",c(2:8))
  
  p_frame <- rbind(p_frame,do.call('c',p))
}

save.image("./reproducibility/mixture_model_with_encoded_features.RData")


colnames(p_frame) <- c("Original",c(2:8))

pdf("./figures/mixture_model_with_encoded_features.pdf",8,4.5)
boxplot(p_frame,xlab="Number of Features", ylab="Accuracy (ARI)", main= "Mixture model with encoded feature information")
dev.off()







