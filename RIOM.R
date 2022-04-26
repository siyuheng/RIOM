#' Function for calculating warning accuracy and sensitivity weights in randomization inference with outcome misclassification
#' @param treat.ind A N-Length vector of the treatment indicators of all N units in the
#'          study, 1 if treated and 0 if not.
#' @param outcome.ind A N-length vector of the measured binary outcomes.
#' @param index A N-length vector: the n-th entry is the index of the stratum of unit n.
#' @param alpha The alpha level of the two-sided test. The default is 0.05.
#' @param null.hypothesis The null hypothesis of interest. Two options: "sharp" = Fisher's sharp null; "weak" = Neyman's weak null.
#' @param type.random The type of the randomization design. Two options: "1" = a type 1 randomization design; "2" = a type 2 randomization design.
#' @param timelimit The limit for the runtime in seconds of the function. The default is 1000 seconds.
#' @param gap The tolerable gap of the upper bound and lower bound of the optimial solution. The default is 0.00001.
#' @return A list that contains: "p-value" - two-sided p-value based on measured outcomes;
#'                               "Difference-in-means Estimate" - the point estimate of the average treatment effect based on measured outcomes; only for weak null.
#'                               "Confidence Interval" - the confidence interval of the average treatment effect based on measured outcomes; only for weak null.
#'                               "Warning Accuracy" - the warning accuracy given the observed data and alpha level.
#'                               "Minimal Alteration Number" - the minimal alteration number given the observed data and alpha level.
#'                               "Sensitivity Weights" - the four sensitivity weights given the observed data and alpha level.
#'                               "Runtime" - the total run time in seconds.
#' @import gurobi, mgcv, slam
#' @export


RIOM<-function(treat.ind, outcome.ind, index, alpha = 0.05, null.hypothesis, type.random, timelimit = 1000,  gap = 0.00001){
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package \"gurobi\" is required for this function.", call. = FALSE)
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package \"mgcv\" is required for this function.", call. = FALSE)
  }
  start_time <- Sys.time()
  I=max(index)
  Z<-rep(0, length(treat.ind))
  Y<-rep(0, length(outcome.ind))
  n<-rep(0, I)
  m<-rep(0, I)
  Q<-rep(0, length(Z))
  Q_sub<-rep(0, length(Z))
  start.point=1
  end.point=sum(as.numeric(index==1))
  for (i in 1:I){
    Z[start.point:end.point]=treat.ind[index==i]
    Y[start.point:end.point]=outcome.ind[index==i]
    m[i]=sum(as.numeric(treat.ind[index==i]==1))
    n[i]=sum(as.numeric(index==i))
    Q[start.point:end.point]=i
    Q_sub[start.point:end.point]=c(1:sum(as.numeric(index==i)))
    if (i<I){
      start.point=end.point+1
      end.point=end.point+sum(as.numeric(index==(i+1)))
    }
  }
  c=qchisq(1-alpha, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE)

  ############# Sharp Type 1 #############################

  if (null.hypothesis=="sharp" & type.random==1){
    T_MH=sum(Z*Y)
    E_null=0
    V_null=0
    for (i in 1:I){
      E_null=E_null+(m[i]/n[i])*sum(Y[Q==i])
      V_null=V_null+(m[i]*sum(Y[Q==i])*(n[i]-sum(Y[Q==i]))*(n[i]-m[i]))/(n[i]^2*((n[i])-1))
    }
    T_MH_normalized<-((T_MH-E_null)^2)/V_null
    pvalue<-pchisq(T_MH_normalized, df=1, lower.tail = FALSE)

    if (pvalue>=alpha){
      my_list=list("p-value" = pvalue, "remark" = "p-value is above the significance level")
      return(my_list)
    }

    K_1<-rep(0, 4*I)
    K_2<-rep(0, 4*I)
    for (i in 1:I){
      K_1[4*i-3]=i
      K_1[4*i-2]=i
      K_1[4*i-1]=i
      K_1[4*i]=i
      K_2[4*i-3]=1
      K_2[4*i-2]=2
      K_2[4*i-1]=3
      K_2[4*i]=4
    }

    q<-rep(0, 4*I)
    for (i in 1:I){
      q[4*i-3]=-1/sum(n)
      q[4*i-2]=1/sum(n)
      q[4*i-1]=-1/sum(n)
      q[4*i]=1/sum(n)
    }

    Q_1<-matrix(0, nrow=4*I, ncol=4*I)
    for (s in 1:(4*I)){
      for (t in 1:(4*I)){
        if (K_1[s]==K_1[t]){
          if ((K_2[s]==1 & K_2[t]==1)|(K_2[s]==1 & K_2[t]==2)|(K_2[s]==2 & K_2[t]==1)|(K_2[s]==2 & K_2[t]==2)){
            Q_1[s,t]=(m[K_1[s]]/n[K_1[s]])^2+c*(m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/(n[K_1[s]]^2*(n[K_1[s]]-1))
          } else if ((K_2[s]==3 & K_2[t]==3)|(K_2[s]==3 & K_2[t]==4)|(K_2[s]==4 & K_2[t]==3)|(K_2[s]==4 & K_2[t]==4)){
            Q_1[s,t]=(1-m[K_1[s]]/n[K_1[s]])^2+c*(m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/(n[K_1[s]]^2*(n[K_1[s]]-1))
          } else {
            Q_1[s,t]=-(m[K_1[s]]/n[K_1[s]])*(1-m[K_1[s]]/n[K_1[s]])+c*(m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/(n[K_1[s]]^2*(n[K_1[s]]-1))
          }
        } else {
          if ((K_2[s]==1 & K_2[t]==1)|(K_2[s]==1 & K_2[t]==2)|(K_2[s]==2 & K_2[t]==1)|(K_2[s]==2 & K_2[t]==2)){
            Q_1[s,t]=(m[K_1[s]]/n[K_1[s]])*(m[K_1[t]]/n[K_1[t]])
          } else if ((K_2[s]==1 & K_2[t]==3)|(K_2[s]==1 & K_2[t]==4)|(K_2[s]==2 & K_2[t]==3)|(K_2[s]==2 & K_2[t]==4)){
            Q_1[s,t]=-(m[K_1[s]]/n[K_1[s]])*(1-m[K_1[t]]/n[K_1[t]])
          } else if ((K_2[s]==3 & K_2[t]==1)|(K_2[s]==3 & K_2[t]==2)|(K_2[s]==4 & K_2[t]==1)|(K_2[s]==4 & K_2[t]==2)){
            Q_1[s,t]=-(1-m[K_1[s]]/n[K_1[s]])*(m[K_1[t]]/n[K_1[t]])
          } else {
            Q_1[s,t]=(1-m[K_1[s]]/n[K_1[s]])*(1-m[K_1[t]]/n[K_1[t]])
          }
        }
      }
    }
    q_1<-rep(0, 4*I)
    for (s in 1:(4*I)){
      q_1[s]=-c*((m[K_1[s]]*n[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/((n[K_1[s]])^2*(n[K_1[s]]-1)))
    }

    ##################################################
    # Initiate a model
    model = list()
    model$modelsense = 'max'
    model$vtype = rep('I', 4*I)
    L<-rep(0, 4*I)
    U<-rep(0, 4*I)
    for (i in 1:I){
      U[4*i-3]=sum((1-Z[Q==i])*(1-Y[Q==i]))
      U[4*i-2]=sum((1-Z[Q==i])*Y[Q==i])
      U[4*i-1]=sum(Z[Q==i]*(1-Y[Q==i]))
      U[4*i]=sum(Z[Q==i]*Y[Q==i])
    }
    model$lb<-L
    model$ub<-U

    ###################################################
    #Objective function
    model$obj = q  #linear term
    model$objcon = sum(1-Y)/sum(n)     #constant term

    ###########################################################
    # Linear constraint
    model$A=matrix(0, nrow=1, ncol=4*I)

    ###########################################################
    # Quadratic constraint
    model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '<')

    #if (verbose) flag = 1
    #else flag = 0

    #https://www.gurobi.com/documentation/9.1/refman/parameters.html#sec:Parameters
    res = gurobi(model, params = list(NonConvex = 2, MIPGapAbs = gap, OutputFlag = 0))

    if (as.numeric(res$status=="OPTIMAL") == 0){
      my_list=list("p-value" = pvalue, "remark" = "Runtime surpassed the prespecified time limit")
      return(my_list)
    }

    WA=res$objval
    AN=(1-res$objval)*sum(n)

    end_time <- Sys.time()

    runtime=total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    #https://www.gurobi.com/documentation/9.1/quickstart_mac/cs_results.html

    ##############Sensitivity Weights########################
    SW=matrix(0, nrow = 2, ncol = 2)
    colnames(SW)<-c("False Positives", "False Negatives")
    row.names(SW)<-c("Treated", "Control")
    SW_00=0
    SW_01=0
    SW_10=0
    SW_11=0

    for (i in 1:I){
      SW_00=SW_00+res$x[4*i-3]
      SW_01=SW_01+U[4*i-2]-res$x[4*i-2]
      SW_10=SW_10+res$x[4*i-1]
      SW_11=SW_11+U[4*i]-res$x[4*i]
    }

    SW[2,2]=SW_00/(SW_00+SW_01+SW_10+SW_11)
    SW[2,1]=SW_01/(SW_00+SW_01+SW_10+SW_11)
    SW[1,2]=SW_10/(SW_00+SW_01+SW_10+SW_11)
    SW[1,1]=SW_11/(SW_00+SW_01+SW_10+SW_11)

    my_list=list("p-value" = pvalue, "Warning Accuracy" = WA, "Minimal Alteration Number" = AN, "Sensitivity Weights" = SW, "Runtime" = total_time)

    return(my_list)

  }

  ################# Sharp Type 2 ################################

  if (null.hypothesis=="sharp" & type.random==2){
    T_MH=sum(Z*Y)
    E_null=0
    V_null=0
    for (i in 1:I){
      E_null=E_null+(m[i]/n[i])*sum(Y[Q==i])
      V_null=V_null+(m[i]*sum(Y[Q==i])*(n[i]-sum(Y[Q==i]))*(n[i]-m[i]))/(n[i]^2*((n[i])-1))
    }
    T_MH_normalized<-((T_MH-E_null)^2)/V_null
    pvalue<-pchisq(T_MH_normalized, df=1, lower.tail = FALSE)

    if (pvalue>=alpha){
      my_list=list("p-value" = pvalue, "remark" = "p-value is above the significance level")
      return(my_list)
    }

    I_info<-matrix(0, nrow = I, ncol = 4) # 1-00, 2-01, 3-10, 4-11
    for (i in 1:I){
      I_info[i,1]=sum((1-Z[Q==i])*(1-Y[Q==i]))
      I_info[i,2]=sum((1-Z[Q==i])*Y[Q==i])
      I_info[i,3]=sum(Z[Q==i]*(1-Y[Q==i]))
      I_info[i,4]=sum(Z[Q==i]*Y[Q==i])
    }
    S_unique<-uniquecombs(I_info,ordered=FALSE) # Number of unique 2 by 2 tables
    S=nrow(S_unique)
    P_s<-rep(0, S)
    for (s in 1:S){
      count=0
      for (i in 1:I){
        if (sum(as.numeric(I_info[i,]==S_unique[s,]))==4){
          count=count+1
        }
      }
      P_s[s]=count
    }
    N_s<-rep(0, S)
    for (s in 1:S){
      N_s[s]=(S_unique[s,1]+1)*(S_unique[s,2]+1)*(S_unique[s,3]+1)*(S_unique[s,4]+1)
    }
    Delta<-matrix(0, nrow = sum(N_s), ncol = 4)
    for (s in 1:S){
      C_00<-seq(from = 0, to = S_unique[s,1], by = 1)
      C_01<-seq(from = 0, to = S_unique[s,2], by = 1)
      C_10<-seq(from = 0, to = S_unique[s,3], by = 1)
      C_11<-seq(from = 0, to = S_unique[s,4], by = 1)
      C_all<-as.matrix(expand.grid(C_00, C_01, C_10, C_11))
      if (s == 1){
        Delta[c(1:N_s[1]),]=C_all
      } else {
        end_1<-sum(N_s[1:(s-1)])+1
        end_2<-sum(N_s[1:s])
        Delta[c(end_1:end_2), ]=C_all
      }
    }
    N_s_index<-rep(0, sum(N_s))
    for (s in 1:S){
      if (s == 1){
        N_s_index[1:N_s[1]]=1
      } else {
        end_1<-sum(N_s[1:(s-1)])+1
        end_2<-sum(N_s[1:s])
        N_s_index[end_1:end_2]=s
      }
    }
    n_s<-rep(0, S)
    m_s<-rep(0, S)
    for (s in 1:S){
      n_s[s]=sum(S_unique[s, ])
      m_s[s]=S_unique[s,3]+S_unique[s,4]
    }

    ##################################################
    # Initiate a model
    model = list()
    model$modelsense = 'max'
    model$vtype = rep('I', sum(N_s))
    L<-rep(0, sum(N_s))
    model$lb<-L

    ###################################################
    #Objective function
    q<-rep(0, sum(N_s))
    for (i in 1:sum(N_s)){
      q[i]=(Delta[i,2]+Delta[i,4]-Delta[i,1]-Delta[i,3])/(sum(n))
    }
    model$obj = q  #linear term
    model$objcon = sum(1-Y)/sum(n)     #constant term

    ###########################################################
    # Linear constraint
    A_linear=matrix(0, nrow=length(N_s), ncol=sum(N_s))
    for (i in 1:length(N_s)){
      if (i == 1){
        A_linear[i, c(1:N_s[1])]=rep(1, N_s[1])
      } else {
        end_1<-sum(N_s[1:(i-1)])+1
        end_2<-sum(N_s[1:i])
        A_linear[i, c(end_1:end_2)]=rep(1, N_s[i])
      }
    }
    model$A=A_linear
    model$rhs=P_s
    model$sense='='

    ###########################################################
    # Quadratic constraint
    Q_1<-matrix(0, nrow = sum(N_s), ncol = sum(N_s))
    Q_1_i<-matrix(0, nrow = sum(N_s), ncol = 1)
    Q_1_j<-matrix(0, nrow = 1, ncol = sum(N_s))
    for (i in 1:sum(N_s)){
      Q_1_i[i,1]=(Delta[i,3]+Delta[i,4])-(m_s[N_s_index[i]]/n_s[N_s_index[i]])*sum(Delta[i, ])
      Q_1_j[1,i]=(Delta[i,3]+Delta[i,4])-(m_s[N_s_index[i]]/n_s[N_s_index[i]])*sum(Delta[i, ])
    }
    Q_1<-Q_1_i%*%Q_1_j
    q_1<-rep(0, sum(N_s))
    for (i in 1:sum(N_s)){
      q_1[i]=-c*(m_s[N_s_index[i]]*sum(Delta[i,])*(n_s[N_s_index[i]]-sum(Delta[i,]))*(n_s[N_s_index[i]]-m_s[N_s_index[i]]))/((n_s[N_s_index[i]])^2*(n_s[N_s_index[i]]-1))
    }

    model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '<')

    #https://www.gurobi.com/documentation/9.1/refman/parameters.html#sec:Parameters
    res = gurobi(model, params = list(NonConvex = 2, MIPGapAbs = gap, OutputFlag = 0, TimeLimit=timelimit))

    if (as.numeric(res$status=="OPTIMAL") == 0){
      my_list=list("p-value" = pvalue, "remark" = "Runtime surpassed the prespecified time limit")
      return(my_list)
    }

    WA=res$objval
    AN=(1-res$objval)*sum(n)

    end_time <- Sys.time()

    runtime=total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))


    ##############Sensitivity Weights########################
    SW=matrix(0, nrow = 2, ncol = 2)
    colnames(SW)<-c("False Positives", "False Negatives")
    row.names(SW)<-c("Treated", "Control")
    SW_00=0
    SW_01=0
    SW_10=0
    SW_11=0

    # S_unique 1-00, 2-01, 3-10, 4-11
    for (i in 1:sum(N_s)){
      SW_00=SW_00+res$x[i]*Delta[i,1]
      SW_01=SW_01+res$x[i]*(S_unique[N_s_index[i], 2]-Delta[i,2])
      SW_10=SW_10+res$x[i]*Delta[i,3]
      SW_11=SW_11+res$x[i]*(S_unique[N_s_index[i], 4]-Delta[i,4])
    }

    SW[2,2]=SW_00/(SW_00+SW_01+SW_10+SW_11)
    SW[2,1]=SW_01/(SW_00+SW_01+SW_10+SW_11)
    SW[1,2]=SW_10/(SW_00+SW_01+SW_10+SW_11)
    SW[1,1]=SW_11/(SW_00+SW_01+SW_10+SW_11)

    my_list=list("p-value" = pvalue, "Warning Accuracy" = WA, "Minimal Alteration Number" = AN, "Sensitivity Weights" = SW, "Runtime" = total_time)

    return(my_list)
  }

  ################ Weak Type 1 #########################

  if (null.hypothesis=="weak" & type.random==1){
    T_Neyman=0
    V_null_1=0
    V_null_2=0
    for (i in 1:I){
      T_Neyman=T_Neyman+(n[i]/sum(n))*((sum(Z[Q==i]*Y[Q==i]))/m[i]-(sum((1-Z[Q==i])*Y[Q==i]))/(n[i]-m[i]))
      V_T=sum(Z[Q==i]*(Y[Q==i]-(1/m[i])*sum(Z[Q==i]*Y[Q==i]))^2)
      V_C=sum((1-Z[Q==i])*(Y[Q==i]-(1/(n[i]-m[i]))*sum((1-Z[Q==i])*Y[Q==i]))^2)
      V_null_1=V_null_1+(n[i]/sum(n))^2*((1/(m[i]*(m[i]-1)))*V_T+(1/((n[i]-m[i])*(n[i]-m[i]-1)))*V_C)
      V_null_2=V_null_2+(n[i]/sum(n))^2*((1/m[i])*var(Y[Q==i & Z==1])+(1/(n[i]-m[i]))*var(Y[Q==i & Z==0]))
    }
    T_Neyman_normalized<-((T_Neyman)^2)/V_null_1
    pvalue<-pchisq(T_Neyman_normalized, df=1, lower.tail = FALSE)
    CI<-c(T_Neyman-qnorm(0.975)*sqrt(V_null_1), T_Neyman+qnorm(0.975)*sqrt(V_null_1))
    names(CI)=c("Left Endpoint", "Right Endpoint")

    if (pvalue>=alpha){
      my_list=list("p-value" = pvalue, "remark" = "p-value is above the significance level")
      return(my_list)
    }

    K_1<-rep(0, 4*I)
    K_2<-rep(0, 4*I)
    for (i in 1:I){
      K_1[4*i-3]=i
      K_1[4*i-2]=i
      K_1[4*i-1]=i
      K_1[4*i]=i
      K_2[4*i-3]=1
      K_2[4*i-2]=2
      K_2[4*i-1]=3
      K_2[4*i]=4
    }

    q<-rep(0, 4*I)
    for (i in 1:I){
      q[4*i-3]=-1/sum(n)
      q[4*i-2]=1/sum(n)
      q[4*i-1]=-1/sum(n)
      q[4*i]=1/sum(n)
    }

    Q_1<-matrix(0, nrow=4*I, ncol=4*I)
    for (s in 1:(4*I)){
      for (t in 1:(4*I)){
        if (K_1[s]==K_1[t]){
          if ((K_2[s]==1 & K_2[t]==1)|(K_2[s]==1 & K_2[t]==2)|(K_2[s]==2 & K_2[t]==1)|(K_2[s]==2 & K_2[t]==2)){
            Q_1[s,t]=(n[K_1[s]]^2)/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])^2)+c*(n[K_1[s]]^2)/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])^2*(n[K_1[s]]-m[K_1[s]]-1))
          } else if ((K_2[s]==3 & K_2[t]==3)|(K_2[s]==3 & K_2[t]==4)|(K_2[s]==4 & K_2[t]==3)|(K_2[s]==4 & K_2[t]==4)){
            Q_1[s,t]=(n[K_1[s]]^2)/((sum(n))^2*(m[K_1[s]])^2)+c*(n[K_1[s]]^2)/((sum(n))^2*(m[K_1[s]])^2*(m[K_1[s]]-1))
          } else {
            Q_1[s,t]=-(n[K_1[s]]^2)/((sum(n))^2*m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))
          }
        } else {
          if ((K_2[s]==1 & K_2[t]==1)|(K_2[s]==1 & K_2[t]==2)|(K_2[s]==2 & K_2[t]==1)|(K_2[s]==2 & K_2[t]==2)){
            Q_1[s,t]=(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])*(n[K_1[t]]-m[K_1[t]]))
          } else if ((K_2[s]==1 & K_2[t]==3)|(K_2[s]==1 & K_2[t]==4)|(K_2[s]==2 & K_2[t]==3)|(K_2[s]==2 & K_2[t]==4)){
            Q_1[s,t]=-(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])*m[K_1[t]])
          } else if ((K_2[s]==3 & K_2[t]==1)|(K_2[s]==3 & K_2[t]==2)|(K_2[s]==4 & K_2[t]==1)|(K_2[s]==4 & K_2[t]==2)){
            Q_1[s,t]=-(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*m[K_1[s]]*(n[K_1[t]]-m[K_1[t]]))
          } else {
            Q_1[s,t]=(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*m[K_1[s]]*m[K_1[t]])
          }
        }
      }
    }
    q_1<-rep(0, 4*I)
    for (s in 1:(4*I)){
      if (K_2[s]==1 | K_2[s]==2){
        q_1[s]=-c*(n[K_1[s]]^2)/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])*(n[K_1[s]]-m[K_1[s]]-1))
      } else {
        q_1[s]=-c*(n[K_1[s]]^2)/((sum(n))^2*m[K_1[s]]*(m[K_1[s]]-1))
      }
    }


    ##################################################
    # Initiate a model
    model = list()
    model$modelsense = 'max'
    model$vtype = rep('I', 4*I)
    L<-rep(0, 4*I)
    U<-rep(0, 4*I)
    for (i in 1:I){
      U[4*i-3]=sum((1-Z[Q==i])*(1-Y[Q==i]))
      U[4*i-2]=sum((1-Z[Q==i])*Y[Q==i])
      U[4*i-1]=sum(Z[Q==i]*(1-Y[Q==i]))
      U[4*i]=sum(Z[Q==i]*Y[Q==i])
    }
    model$lb<-L
    model$ub<-U

    ###################################################
    #Objective function
    model$obj = q  #linear term
    model$objcon = sum(1-Y)/sum(n)     #constant term

    ###########################################################
    # Linear constraint
    model$A=matrix(0, nrow=1, ncol=4*I)

    ###########################################################
    # Quadratic constraint
    model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '<')

    #if (verbose) flag = 1
    #else flag = 0

    #https://www.gurobi.com/documentation/9.1/refman/parameters.html#sec:Parameters
    res = gurobi(model, params = list(NonConvex = 2, MIPGapAbs = gap, OutputFlag = 0, TimeLimit = timelimit))

    if (as.numeric(res$status=="OPTIMAL") == 0){
      my_list=list("p-value" = pvalue, "remark" = "Runtime surpassed the prespecified time limit")
      return(my_list)
    }

    WA=res$objval
    AN=(1-res$objval)*sum(n)

    end_time <- Sys.time()

    runtime=total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    #https://www.gurobi.com/documentation/9.1/quickstart_mac/cs_results.html


    ##############Sensitivity Weights########################
    SW=matrix(0, nrow = 2, ncol = 2)
    colnames(SW)<-c("False Positives", "False Negatives")
    row.names(SW)<-c("Treated", "Control")
    SW_00=0
    SW_01=0
    SW_10=0
    SW_11=0

    for (i in 1:I){
      SW_00=SW_00+res$x[4*i-3]
      SW_01=SW_01+U[4*i-2]-res$x[4*i-2]
      SW_10=SW_10+res$x[4*i-1]
      SW_11=SW_11+U[4*i]-res$x[4*i]
    }

    SW[2,2]=SW_00/(SW_00+SW_01+SW_10+SW_11)
    SW[2,1]=SW_01/(SW_00+SW_01+SW_10+SW_11)
    SW[1,2]=SW_10/(SW_00+SW_01+SW_10+SW_11)
    SW[1,1]=SW_11/(SW_00+SW_01+SW_10+SW_11)

    my_list=list("p-value" = pvalue, 'Difference-in-means Estimate' = T_Neyman, 'Confidence Interval' = CI, "Warning Accuracy" = WA, "Minimal Alteration Number" = AN, "Sensitivity Weights" = SW, "Runtime" = total_time)

    return(my_list)

  }

  if (null.hypothesis=="weak" & type.random==2){
    T_Neyman=0
    V_null_1=0
    V_null_2=0
    for (i in 1:I){
      T_Neyman=T_Neyman+(n[i]/sum(n))*((sum(Z[Q==i]*Y[Q==i]))/m[i]-(sum((1-Z[Q==i])*Y[Q==i]))/(n[i]-m[i]))
      V_T=sum(Z[Q==i]*(Y[Q==i]-(1/m[i])*sum(Z[Q==i]*Y[Q==i]))^2)
      V_C=sum((1-Z[Q==i])*(Y[Q==i]-(1/(n[i]-m[i]))*sum((1-Z[Q==i])*Y[Q==i]))^2)
      V_null_1=V_null_1+(n[i]/sum(n))^2*((1/(m[i]*(m[i]-1)))*V_T+(1/((n[i]-m[i])*(n[i]-m[i]-1)))*V_C)
      V_null_2=V_null_2+(n[i]/sum(n))^2*((1/m[i])*var(Y[Q==i & Z==1])+(1/(n[i]-m[i]))*var(Y[Q==i & Z==0]))
    }
    T_Neyman_normalized<-(T_Neyman^2)/V_null_1
    pvalue<-pchisq(T_Neyman_normalized, df=1, lower.tail = FALSE)
    CI<-c(T_Neyman-qnorm(0.975)*sqrt(V_null_1), T_Neyman+qnorm(0.975)*sqrt(V_null_1))
    names(CI)=c("Left Endpoint", "Right Endpoint")

    if (pvalue>=alpha){
      my_list=list("p-value" = pvalue, "remark" = "p-value is above the significance level")
      return(my_list)
    }

    I_info<-matrix(0, nrow = I, ncol = 4) # 1-00, 2-01, 3-10, 4-11
    for (i in 1:I){
      I_info[i,1]=sum((1-Z[Q==i])*(1-Y[Q==i]))
      I_info[i,2]=sum((1-Z[Q==i])*Y[Q==i])
      I_info[i,3]=sum(Z[Q==i]*(1-Y[Q==i]))
      I_info[i,4]=sum(Z[Q==i]*Y[Q==i])
    }
    S_unique<-uniquecombs(I_info,ordered=FALSE) # Number of unique 2 by 2 tables
    S=nrow(S_unique)
    P_s<-rep(0, S)
    for (s in 1:S){
      count=0
      for (i in 1:I){
        if (sum(as.numeric(I_info[i,]==S_unique[s,]))==4){
          count=count+1
        }
      }
      P_s[s]=count
    }
    N_s<-rep(0, S)
    for (s in 1:S){
      N_s[s]=(S_unique[s,1]+1)*(S_unique[s,2]+1)*(S_unique[s,3]+1)*(S_unique[s,4]+1)
    }
    Delta<-matrix(0, nrow = sum(N_s), ncol = 4)
    for (s in 1:S){
      C_00<-seq(from = 0, to = S_unique[s,1], by = 1)
      C_01<-seq(from = 0, to = S_unique[s,2], by = 1)
      C_10<-seq(from = 0, to = S_unique[s,3], by = 1)
      C_11<-seq(from = 0, to = S_unique[s,4], by = 1)
      C_all<-as.matrix(expand.grid(C_00, C_01, C_10, C_11))
      if (s == 1){
        Delta[c(1:N_s[1]),]=C_all
      } else {
        end_1<-sum(N_s[1:(s-1)])+1
        end_2<-sum(N_s[1:s])
        Delta[c(end_1:end_2), ]=C_all
      }
    }
    N_s_index<-rep(0, sum(N_s))
    for (s in 1:S){
      if (s == 1){
        N_s_index[1:N_s[1]]=1
      } else {
        end_1<-sum(N_s[1:(s-1)])+1
        end_2<-sum(N_s[1:s])
        N_s_index[end_1:end_2]=s
      }
    }
    n_s<-rep(0, S)
    m_s<-rep(0, S)
    for (s in 1:S){
      n_s[s]=sum(S_unique[s, ])
      m_s[s]=S_unique[s,3]+S_unique[s,4]
    }

    ##################################################
    # Initiate a model
    model = list()
    model$modelsense = 'max'
    model$vtype = rep('I', sum(N_s))
    L<-rep(0, sum(N_s))
    model$lb<-L

    ###################################################
    #Objective function
    q<-rep(0, sum(N_s))
    for (i in 1:sum(N_s)){
      q[i]=(Delta[i,2]+Delta[i,4]-Delta[i,1]-Delta[i,3])/(sum(n))
    }
    model$obj = q  #linear term
    model$objcon = sum(1-Y)/sum(n)     #constant term

    ###########################################################
    # Linear constraint
    A_linear=matrix(0, nrow=length(N_s), ncol=sum(N_s))
    for (i in 1:length(N_s)){
      if (i == 1){
        A_linear[i, c(1:N_s[1])]=rep(1, N_s[1])
      } else {
        end_1<-sum(N_s[1:(i-1)])+1
        end_2<-sum(N_s[1:i])
        A_linear[i, c(end_1:end_2)]=rep(1, N_s[i])
      }
    }
    model$A=A_linear
    model$rhs=P_s
    model$sense='='

    ###########################################################
    # Quadratic constraint
    Q_1<-matrix(0, nrow = sum(N_s), ncol = sum(N_s))
    Q_1_i<-matrix(0, nrow = sum(N_s), ncol = 1)
    Q_1_j<-matrix(0, nrow = 1, ncol = sum(N_s))
    for (i in 1:sum(N_s)){
      Q_1_i[i,1]=(n_s[N_s_index[i]]*(Delta[i, 3]+Delta[i,4]))/(sum(n)*m_s[N_s_index[i]])-(n_s[N_s_index[i]]*(Delta[i, 1]+Delta[i,2]))/(sum(n)*(n_s[N_s_index[i]]-m_s[N_s_index[i]]))
      Q_1_j[1,i]=(n_s[N_s_index[i]]*(Delta[i, 3]+Delta[i,4]))/(sum(n)*m_s[N_s_index[i]])-(n_s[N_s_index[i]]*(Delta[i, 1]+Delta[i,2]))/(sum(n)*(n_s[N_s_index[i]]-m_s[N_s_index[i]]))
    }
    Q_1<-Q_1_i%*%Q_1_j
    q_1<-rep(0, sum(N_s))
    for (i in 1:sum(N_s)){
      q_1_1=(Delta[i,3]+Delta[i,4])/(m_s[N_s_index[i]]*(m_s[N_s_index[i]]-1))
      q_1_2=((Delta[i,3]+Delta[i,4])^2)/(m_s[N_s_index[i]]^2*(m_s[N_s_index[i]]-1))
      q_1_3=(Delta[i,1]+Delta[i,2])/((n_s[N_s_index[i]]-m_s[N_s_index[i]])*(n_s[N_s_index[i]]-m_s[N_s_index[i]]-1))
      q_1_4=((Delta[i,1]+Delta[i,2])^2)/(((n_s[N_s_index[i]]-m_s[N_s_index[i]])^2)*(n_s[N_s_index[i]]-m_s[N_s_index[i]]-1))
      q_1[i]=-c*((n_s[N_s_index[i]]/sum(n))^2)*(q_1_1-q_1_2+q_1_3-q_1_4)
    }

    model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '<')

    #if (verbose) flag = 1
    #else flag = 0
    #https://www.gurobi.com/documentation/9.1/refman/parameters.html#sec:Parameters
    res = gurobi(model, params = list(NonConvex = 2, MIPGapAbs = gap, OutputFlag = 0, TimeLimit=timelimit))

    if (as.numeric(res$status=="OPTIMAL") == 0){
      my_list=list("p-value" = pvalue, "remark" = "Runtime surpassed the prespecified time limit")
      return(my_list)
    }

    WA=res$objval
    AN=(1-res$objval)*sum(n)

    end_time <- Sys.time()

    runtime=total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    ##############Sensitivity Weights########################
    SW=matrix(0, nrow = 2, ncol = 2)
    colnames(SW)<-c("False Positives", "False Negatives")
    row.names(SW)<-c("Treated", "Control")
    SW_00=0
    SW_01=0
    SW_10=0
    SW_11=0

    # S_unique 1-00, 2-01, 3-10, 4-11
    for (i in 1:sum(N_s)){
      SW_00=SW_00+res$x[i]*Delta[i,1]
      SW_01=SW_01+res$x[i]*(S_unique[N_s_index[i], 2]-Delta[i,2])
      SW_10=SW_10+res$x[i]*Delta[i,3]
      SW_11=SW_11+res$x[i]*(S_unique[N_s_index[i], 4]-Delta[i,4])
    }

    SW[2,2]=SW_00/(SW_00+SW_01+SW_10+SW_11)
    SW[2,1]=SW_01/(SW_00+SW_01+SW_10+SW_11)
    SW[1,2]=SW_10/(SW_00+SW_01+SW_10+SW_11)
    SW[1,1]=SW_11/(SW_00+SW_01+SW_10+SW_11)

    my_list=list("p-value" = pvalue, 'Difference-in-means Estimate' = T_Neyman, 'Confidence Interval' = CI, "Warning Accuracy" = WA, "Minimal Alteration Number" = AN, "Sensitivity Weights" = SW, "Runtime" = total_time)

    return(my_list)

  }
}


