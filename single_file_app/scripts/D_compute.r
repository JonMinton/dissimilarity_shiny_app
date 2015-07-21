     Dissimilarity.compute <- function(minority, total, n.boot=1000)
     {
     #### Compute the dissimilarity index
     proportion <- minority / total
     proportion.overall <- sum(minority) / sum(total)
     Total <- sum(total)
     D <- sum(total * abs(proportion - proportion.overall)) / (2 * Total * proportion.overall * (1-proportion.overall))    
     
     #### Compute bootstrapped estimates of uncertainty
     D.boot <- rep(NA, n.boot)
     n <- length(total)
          for(j in 1:n.boot)
          {
          ## Generate a bootstrap sample
          samp <- sample(x=1:n, size=n, replace=TRUE)
          minority.boot <- minority[samp]
          total.boot <- total[samp]
          
          ## Compute D based on this sample
          proportion.boot <- minority.boot / total.boot
          proportion.overall.boot <- sum(minority.boot) / sum(total.boot)
          Total.boot <- sum(total.boot)
          D.boot[j] <- sum(total.boot * abs(proportion.boot - proportion.overall.boot)) / (2 * Total.boot * proportion.overall.boot * (1-proportion.overall.boot))    
          }
     D.estimate <- c(D, quantile(D.boot, c(0.025, 0.5, 0.975), na.rm=T))
     results <- list(D.estimate=D.estimate, D.boot=D.boot)
     return(results)
     }
