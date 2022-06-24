
source("./GMMAT_helper_functions.R")
# devtools::install_github("AskExplain/gcode_R@beta_test_v2022.2")

library(quadprog)
library(GMMAT)
library(gcode)

# SNPs
p = 300

# Samples
n = 1000

# Mean Allele Frequency
f = 0.5



time_run_yes_200_encode <- c()
theta_param_yes_200_encode <- c()

time_run_no_encode <- c()
theta_param_no_encode <- c()

time_run_HE_encode <- c()
theta_param_HE_encode <- c()


# Heritability at 0.5
for(h2 in c(0.5)){
  
  # Number of permutations
  for (i in 1:30){
    print(i)
    set.seed(i+10*h2)
    
    # Create simulated data per permutation
    X.s = (replicate(p, rbinom(n, size = 2, p = f) )) #use d genotypes
    #tau^2 = var(lambda.s) = h2 / p
    lambda.s = array(rnorm(p, 0, sqrt( h2 / p)),c(p)) #d effects
    #generate phenotype: SNP effects + random noise with var=1-h2
    y = ( X.s %*% lambda.s + rnorm(n, 0, sqrt(1-h2) )) #scaling makes mean(y)=0 as our LMM ignores int#test individual SNPs and make a QQ-plot
    s1 <- Sys.time()
    GRM <- X.s%*%t(X.s)/p
    s2 <- Sys.time()
    scaleX.s <- scaleW(X.s)
    scaleGRM <- scaleX.s%*%t(scaleX.s)/p
    s3 <- Sys.time()
    
    # Run original mixed model
    main_cpu_time_no_encode <- system.time(glmmkin_no_encode <- david_glmmkin.ai(fit0 = glm((y)~1), kins = list(GRM), k = 0, verbose = F, encode = F))

    theta_param_no_encode <- rbind(theta_param_no_encode,
                                   data.frame(iter=i,h2=h2,sigma_e=glmmkin_no_encode$theta[1],sigma_g=glmmkin_no_encode$theta[2])
    )

    time_run_no_encode <- rbind(time_run_no_encode,c(i,h2,main_cpu_time_no_encode,s2-s1))

    
    
    
    
    # Run original mixed model
    main_cpu_time_HE_encode <- system.time(glmmkin_HE_encode <- MASS::ginv(cbind(rbind(c(n*(n-1)/2),sum(scaleGRM)),rbind(sum(scaleGRM),sum(scaleGRM*t(scaleGRM)))))%*%rbind(sum(scale(y)%*%t(scale(y))),sum(scaleGRM*(scale(y)%*%t(scale(y))))))

    theta_param_HE_encode <- rbind(theta_param_HE_encode,
                                   data.frame(iter=i,h2=h2,sigma_e=1-(glmmkin_HE_encode[2]),sigma_g=(glmmkin_HE_encode[2]))
    )
    
    time_run_HE_encode <- rbind(time_run_HE_encode,c(i,h2,main_cpu_time_HE_encode,s3-s1))
    
    
    
    
    
    
    # Run encoded mixed model
    run_encoded_lmm <- function(T){

      k <- 50
      config <- gcode::extract_config(F)
      config$init <- list(alpha="runif",beta="runif")
      config$i_dim <- k
      config$j_dim <- k
      config$verbose <- F
      config$tol <- 1e-3
      config$max_iter <- 1000
      config$n.cores <- 8
      config$learn_rate <- 0.8
      config$batch_size <- 25

      join <- gcode::extract_join_framework(F)
      join$complete <- lapply(join$complete,function(X){c(1)})
      join$covariance <- c(1)

      gcode.model <- gcode::gcode(data_list = list(GRM), config = config, join = join)
      sample_encode <- gcode.model$main.parameters$alpha_sample[[1]]

      glmmkin_yes_encode <- david_glmmkin.ai(fit0 = glm((y)~1), kins = list(GRM), k = k, verbose = F, encode = T, sample_encode = sample_encode)

    }

    main_cpu_time_yes_200_encode <- system.time(

      glmmkin_yes_200_encode <- run_encoded_lmm(T)

    )

    theta_param_yes_200_encode <- rbind(theta_param_yes_200_encode,
                                        data.frame(iter=i,h2=h2,sigma_e=glmmkin_yes_200_encode$theta[1],sigma_g=glmmkin_yes_200_encode$theta[2])
    )

    time_run_yes_200_encode <- rbind(time_run_yes_200_encode,c(i,h2,main_cpu_time_yes_200_encode,s2-s1))
    
    
    
    
    
    
    
    
    
    
    print(
      list(
        
        mean_yes_200 = colMeans(theta_param_yes_200_encode),
        mean_no = colMeans(theta_param_no_encode),
        mean_HE = colMeans(theta_param_HE_encode),
        
        var_yes_200 = apply(theta_param_yes_200_encode,2,var),
        var_no = apply(theta_param_no_encode,2,var),
        var_HE = apply(theta_param_HE_encode,2,var),
        
        time_yes_200 = colMeans(time_run_yes_200_encode),
        time_no = colMeans(time_run_no_encode),
        time_HE = colMeans(time_run_HE_encode)
        
      )
    )
    
  }
  
}

save.image("./reproducibility/encode_original_mixed_model.RData")


pdf("./figures/encoded_vs_original_mixed_model.pdf",width = 15, height = 5)
par(mfcol=c(1,2))
boxplot(data.frame(Encode = theta_param_yes_200_encode[,4], Original = theta_param_no_encode[,4], HE = theta_param_HE_encode[,4]),main="Heritability (h2)",ylim=c(0,1), outline = F)
boxplot(data.frame(Encode = time_run_yes_200_encode[,5], Original = time_run_no_encode[,5], HE = time_run_HE_encode[,5]),main="Runtime (seconds)", outline = F)
dev.off()

