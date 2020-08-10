MAPCO_MLE <- function(by,bx,se_y,se_x,int_beta,int_sigma_alpha,H0="F",b0=0,tol=1e-8,n_iter=3000){

  by = sign(bx)*by
  bx = abs(bx)
  var_y = se_y^2
  var_x = se_x^2
  J=length(by)

  beta = int_beta
  mu_eta = rbind(mean(bx),0)
  int_rho = 0
  Sigma_eta = matrix(c(sd(bx)^2,
                       int_rho*sd(bx)*int_sigma_alpha,
                       int_rho*sd(bx)*int_sigma_alpha,
                       int_sigma_alpha^2),2)
  b = rbind(beta,1)
  c = rbind(1,0)
  iv.Sigma_eta = solve(Sigma_eta)

  maxlnL = diff_lnL = 1000
  iter=0

  while (diff_lnL>tol){

    # posterior distribution of eta_j

    Sigma.p = eta.p = numeric()

    iv_sigma.p = ( matrix(1,J,1)%*%as.vector(iv.Sigma_eta) +
                     as.matrix(1/var_y)%*%as.vector(b%*%t(b)) + as.matrix(1/var_x)%*%as.vector(c%*%t(c)) )
    det_sigma.p = iv_sigma.p[,1]*iv_sigma.p[,4]-(iv_sigma.p[,2])^2
    Sigma.p = (1/det_sigma.p)* cbind(iv_sigma.p[,4],-iv_sigma.p[,2],-iv_sigma.p[,3],iv_sigma.p[,1])

    eta.p_right = t( rep(1,J)%*%(t(mu_eta)%*%iv.Sigma_eta) +
                       (by/var_y)%*%t(b) + (bx/var_x)%*%t(c) )

    eta.p = rbind(colSums(t(Sigma.p[,1:2])*eta.p_right),
                  colSums(t(Sigma.p[,3:4])*eta.p_right))

    # M-step

    if (H0=="F"){
      nub = sum( (-by*eta.p[1,] + eta.p[2,]*eta.p[1,] + Sigma.p[,2]) /var_y )
      deb = sum( ( eta.p[1,]^2 + Sigma.p[,1] )/(-var_y) )
      beta = nub/deb
      b = rbind(beta,1)
    }
    if (H0=="T"){
      beta = b0
      b = rbind(beta,1)
    }

    mu = sum( eta.p[1,])/J
    mu = as.vector(mu)
    mu_eta = rbind(mu,0)

    Sigma_eta = (    rbind(eta.p[1,]-mu_eta[1],eta.p[2,]-mu_eta[2])
                     %*%cbind(eta.p[1,]-mu_eta[1],eta.p[2,]-mu_eta[2])
                     + matrix(rep(1,J)%*%Sigma.p, 2,2) )/J

    iv.Sigma_eta = solve(Sigma_eta)

    # log-likelihood

    D = matrix(NA,2,2)
    D[1,1] = t(b)%*%Sigma_eta%*%b
    D[1,2]=D[2,1] = t(b)%*%Sigma_eta%*%c
    D[2,2] = t(c)%*%Sigma_eta%*%c

    Di = matrix(1,J,1) %*% as.vector(D) + cbind(var_y,0,0,var_x)
    det_Di = Di[,1]*Di[,4]-Di[,2]^2
    iv.Di = (1/det_Di)* cbind(Di[,4],-Di[,2],-Di[,3],Di[,1])

    M = c(t(b)%*%mu_eta,t(c)%*%mu_eta)
    Mi = cbind(by-M[1], bx-M[2])

    maxlnL_ = sum(log( (det_Di)^(-0.5) * exp(-0.5 *(Mi[,1]^2*iv.Di[,1]
                                                    + 2*Mi[,1]*Mi[,2]*iv.Di[,2]
                                                    + Mi[,2]^2*iv.Di[,4]) )
    ))

    diff_lnL = abs(maxlnL_-maxlnL)
    maxlnL = maxlnL_

    iter = iter+1
    if(iter > n_iter){
      break
    }
  }

  rho_hat = Sigma_eta[1,2]/sqrt(Sigma_eta[1,1])/sqrt(Sigma_eta[2,2])
  sigma_alpha = sqrt(Sigma_eta[2,2])

  return (cbind(beta, sigma_alpha, rho_hat, iter, maxlnL))
}


MAPCO <- function(by,bx,se_y,se_x,int_beta=0,int_sigma_alpha=NA,b0=0,tol=1e-8,n_iter=3000){

  if(is.na(int_sigma_alpha)){int_sigma_alpha = sd(bx)
  }else if (int_sigma_alpha==0){int_sigma_alpha = sd(bx)}

  MAPCO_H0 = MAPCO_MLE (by,bx,se_y,se_x,int_beta,int_sigma_alpha,H0="T",b0,tol,n_iter)
  MAPCO_H1 = MAPCO_MLE (by,bx,se_y,se_x,int_beta,int_sigma_alpha,H0="F",b0,tol,n_iter)
  lrt = 2*(MAPCO_H1[5] - MAPCO_H0[5])
  pval = pchisq(lrt,1,lower.tail = FALSE)

  beta = MAPCO_H1[1]
  sigma_alpha = MAPCO_H1[2]
  rho = MAPCO_H1[3]
  iter = MAPCO_H1[4]
  se_beta = abs(beta/sqrt(lrt))

  dichotomy<-function(f,a,b,c,tol_ci=1e-5){
    if(f(a,c)*f(b,c)>0)
      return(fail="find root is fail!")
    else{
      repeat{
        if(abs(b-a)<tol_ci) break;
        x<-(a+b)/2
        if(f(a,c)*f(x,c)<0) b<-x else a<-x
      }
      return(root=(a+b)/2)
    }
  }

  f <- function(b0,MAPCO_H1){
    MAPCO_H0 = MAPCO_MLE (by,bx,se_y,se_x,int_beta,int_sigma_alpha,H0="T",b0,tol,n_iter)
    lrt = 2*(MAPCO_H1[5] - MAPCO_H0[5])
    pval = pchisq(lrt,1,lower.tail = FALSE)
    return(pval-0.05)
  }

  ul = dichotomy(f,beta+se_beta*0.5, beta+se_beta*5, MAPCO_H1)
  ll = dichotomy(f,beta-se_beta*0.5, beta-se_beta*5, MAPCO_H1)
  CI = c(ll,ul)

  return(list(beta=beta, sigma_alpha=sigma_alpha, rho=rho,iter=iter, pval=pval,CI=CI))

}



