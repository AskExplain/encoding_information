
source("~/Documents/main_files/AskExplain/lmform/Rscripts/lmform_heritability/0_explore/helper_functions.R")

# devtools::document("~/Documents/main_files/AskExplain/generative_encoder/package/beta_test_v2022.2/gcode/")
# devtools::install_local("~/Documents/main_files/AskExplain/generative_encoder/package/beta_test_v2022.2/gcode/",force=T)
# library(coxme)
# library(GMMAT)
# library(gcode)

p = 100
n = 300
f = 0.5

time_run_yes_encode <- c()
theta_param_yes_encode <- c()

time_run_no_encode <- c()
theta_param_no_encode <- c()

for(h2 in c(0.5)){
  
  for (i in 1:1000){
    print(i)
    set.seed(i+10*h2)
    X.s = (replicate(p, rbinom(n, size = 2, p = f) )) #use d genotypes
    #tau^2 = var(lambda.s) = h2 / p
    lambda.s = array(rnorm(p, 0, sqrt( h2 / p)),c(p)) #d effects
    #generate phenotype: SNP effects + random noise with var=1-h2
    y = ( X.s %*% lambda.s + rnorm(n, 0, sqrt(1-h2) )) #scaling makes mean(y)=0 as our LMM ignores int#test individual SNPs and make a QQ-plot
    GRM <- X.s%*%t(X.s)/p
    
    main_cpu_time_no_encode <- system.time(glmmkin_no_encode <- david_glmmkin.ai(fit0 = glm(y~1), kins = list(GRM), k = 0, verbose = F, encode = F,tol=0.01, maxiter = 30))
    
    theta_param_no_encode <- rbind(theta_param_no_encode,
                                   data.frame(iter=i,h2=h2,sigma_e=glmmkin_no_encode$theta[1],sigma_g=glmmkin_no_encode$theta[2])
    )
    
    time_run_no_encode <- rbind(time_run_no_encode,c(i,h2,main_cpu_time_no_encode))
    
    
    run_encoded_lmm <- function(T){
      
      k <- 200
      config <- gcode::extract_config(F)
      config$init <- list(alpha="rnorm",beta="rnorm")
      config$i_dim <- k
      config$j_dim <- k
      config$verbose <- F
      config$tol <- 1e-5
      config$max_iter <- 100
      config$n.cores <- 8
      config$learn_rate <- 0.8
      config$batch_size <- 100
      
      join <- gcode::extract_join_framework(F)
      join$complete <- lapply(join$complete,function(X){c(1)})
      join$covariance <- c(0)
      
      gcode.model <- gcode::gcode(data_list = list(cbind(y,X.s)), config = config, join = join)
      sample_encode <- gcode.model$main.parameters$alpha_sample[[1]]
      
      glmmkin_yes_encode <- david_glmmkin.ai(fit0 = glm(y~1), kins = list(GRM), k = k, verbose = F, encode = T, sample_encode = sample_encode, tol = 0.01, maxiter = 30)
      
    }
    
    main_cpu_time_yes_encode <- system.time(
      
      glmmkin_yes_encode <- run_encoded_lmm(T)
      
    )
    
    theta_param_yes_encode <- rbind(theta_param_yes_encode,
                                    data.frame(iter=i,h2=h2,sigma_e=glmmkin_yes_encode$theta[1],sigma_g=glmmkin_yes_encode$theta[2])
    )
    
    time_run_yes_encode <- rbind(time_run_yes_encode,c(i,h2,main_cpu_time_yes_encode))
    
    
    print(
      list(
        
        mean_yes = colMeans(theta_param_yes_encode),
        mean_no = colMeans(theta_param_no_encode),
        
        var_yes = apply(theta_param_yes_encode,2,var),
        var_no = apply(theta_param_no_encode,2,var),
        
        time_yes = colMeans(time_run_yes_encode),
        time_no = colMeans(time_run_no_encode)
        
      )
    )
    
  }
  
}

save.image("~/Documents/main_files/AskExplain/gcta_summaries/data/workflow/lmm/lmm_fast_accuracy_rerun.RData")


t.test(x = (0.5-theta_param_yes_encode$sigma_g)^2, y = (0.5-theta_param_no_encode$sigma_g)^2, paired = T, var.equal = T)
t.test(x = (0.5-theta_param_yes_encode$sigma_e)^2, y = (0.5-theta_param_no_encode$sigma_e)^2, paired = T, var.equal = T)


