#' Single sample binary helper functions
#'
#' @param DesignTable the design table
#' @param TargetRate target rate
#'
#' @return These functions are used by other functions
#' @export
#'
# ocTable_multi <- function(DesignTable,TargetRate){
#   # browser()
#   ocProbMat <- probmat.func(DesignTable,TargetRate)
#   MAT <- as_tibble(ocProbMat)
#   colnames(MAT) <- paste('Analysis',1:ncol(MAT))
#
#   fgo <- DesignTable$MinGo[nrow(DesignTable)]
#   fnogo <- DesignTable$MaxNoGo[nrow(DesignTable)]
#   N <- DesignTable$n[nrow(DesignTable)]
#
#   TBL <- cbind(X=0:(nrow(MAT)-1),MAT) %>%
#     gather('Analysis','Prob',-X) %>%
#     left_join(DesignTable %>%
#                 mutate(Analysis = paste('Analysis',1:n())),
#               by='Analysis') %>%
#     group_by(Analysis) %>%
#     summarise(`First Go` = sum(Prob[X >= MinGo]),
#               `First No Go` = sum(Prob[X <= MaxNoGo]),
#               `Grey Area` = sum(Prob[X > MaxNoGo & X < MinGo]),
#               `Discordant Go` =`First Go` * Discord_Go(MinGo[1],n[1],fgo,N,TargetRate,Weights = pull(MAT,Analysis[1])[MinGo[1]:n[1]+1]),
#               `Discordant No Go` =`First No Go` * Discord_NoGo(MaxNoGo[1],n[1],fnogo,N,TargetRate,Weights = pull(MAT,Analysis[1])[0:MaxNoGo[1]+1]),
#               `Discordant Grey`= `Grey Area` * Discord_Grey(MaxNoGo[1]+1,MinGo[1]-1,n[1],fgo,fnogo,N,TargetRate,Weights = pull(MAT,Analysis[1])[(MaxNoGo[1]+1):(MinGo[1]-1)+1])
#     )
#
#   TBL %>%
#     rbind(
#       tibble(
#         Analysis = 'Final-solo',
#         `First Go` = 1-pbinom(fgo-1,N,TargetRate),
#         `First No Go` = pbinom(fnogo,N,TargetRate),
#         `Grey Area` = 1-`First Go`-`First No Go`,
#         `Discordant Go` = NA,
#         `Discordant No Go` = NA,
#         `Discordant Grey` = NA
#       )
#     ) %>%
#     rbind(
#       TBL %>%
#         gather('Var','Val',-Analysis) %>%
#         group_by(Var) %>%
#         summarise(Sum = sum(Val)) %>%
#         spread(Var,Sum) %>%
#         mutate(
#           `Grey Area`=TBL$`Grey Area`[nrow(TBL)],
#           `Discordant Grey` =TBL$`Discordant Grey`[nrow(TBL)-1],
#           Analysis='Overall')
#     ) %>%
#     mutate(targetRate=TargetRate)
# }


# GC Update: 8/3/21
#' ocTable_multi
#'
#' @param DesignTable DesignTable
#' @param TargetRate TargetRate
#'
#' @return These functions are used by other functions
#' @export
#'
ocTable_multi <- function(DesignTable, TargetRate){
  # browser()
  ocProbMat <- probmat.func(DesignTable,TargetRate)
  MAT <- as_tibble(ocProbMat)
  colnames(MAT) <- paste('Analysis',1:ncol(MAT))

  fgo <- DesignTable$MinGo[nrow(DesignTable)]
  fnogo <- DesignTable$MaxNoGo[nrow(DesignTable)]
  N <- DesignTable$n[nrow(DesignTable)]

  TBL <- cbind(X=0:(nrow(MAT)-1),MAT) %>%
    gather('Analysis','Prob',-X) %>%
    left_join(DesignTable %>%
                mutate(Analysis = paste('Analysis',1:n())),
              by='Analysis') %>%
    group_by(Analysis) %>%
    summarise(`First Go` = sum(Prob[X >= MinGo]),
              `First No Go` = sum(Prob[X <= MaxNoGo]),
              `Grey Area` = sum(Prob[X > MaxNoGo & X < MinGo]),
              `Discordant Go` =`First Go` * Discord_Go(MinGo[1],n[1],fgo,N,TargetRate,Weights = pull(MAT,Analysis[1])[MinGo[1]:n[1]+1]),
              `Discordant No Go` =`First No Go` * Discord_NoGo(MaxNoGo[1],n[1],fnogo,N,TargetRate,Weights = pull(MAT,Analysis[1])[0:MaxNoGo[1]+1]),
              `Discordant Grey`= ifelse(class(try(`Grey Area` * Discord_Grey(MaxNoGo[1]+1,MinGo[1]-1,n[1],fgo,fnogo,N,TargetRate,Weights = pull(MAT,Analysis[1])[(MaxNoGo[1]+1):(MinGo[1]-1)+1])))=="try-error", NA,
                                        `Grey Area` * Discord_Grey(MaxNoGo[1]+1,MinGo[1]-1,n[1],fgo,fnogo,N,TargetRate,Weights = pull(MAT,Analysis[1])[(MaxNoGo[1]+1):(MinGo[1]-1)+1]))
    )

  TBL %>%
    rbind(
      tibble(
        Analysis = 'Final-solo',
        `First Go` = 1-pbinom(fgo-1,N,TargetRate),
        `First No Go` = pbinom(fnogo,N,TargetRate),
        `Grey Area` = 1-`First Go`-`First No Go`,
        `Discordant Go` = NA,
        `Discordant No Go` = NA,
        `Discordant Grey` = NA
      )
    ) %>%
    rbind(
      TBL %>%
        gather('Var','Val',-Analysis) %>%
        group_by(Var) %>%
        summarise(Sum = sum(Val)) %>%
        spread(Var,Sum) %>%
        mutate(
          `Grey Area`=TBL$`Grey Area`[nrow(TBL)],
          `Discordant Grey` =TBL$`Discordant Grey`[nrow(TBL)-1],
          Analysis='Overall')
    ) %>%
    mutate(targetRate=TargetRate)
}



###########################
#' probmat.func
#'
#' @param DesignTable DesignTable
#' @param TargetRate TargetRate
#'
#' @return These functions are used by other functions
#' @export
#'
probmat.func <- function(DesignTable,TargetRate){
  R <- cbind(
    pmax(DesignTable$MaxNoGo +1,0),
    pmin(DesignTable$MinGo -1,DesignTable$n)
  )
  M <- c(DesignTable$n[1],diff(DesignTable$n))

  Intr_func(M,R,TargetRate)
}
#################################
#' Intr_func
#'
#' @param M M
#' @param R R
#' @param ORR ORR
#'
#' @return These functions are used by other functions
#' @export
#'
Intr_func <- function(M,R,ORR){
  Sizes <- unique(M)
  N <- sum(M)
  probs <- lapply(Sizes,function(n) dbinom(0:n,n,ORR))

  MAT <- matrix(0,nrow=N+1,ncol=length(M))

  MAT[1:(M[1]+1),1] <- probs[[which(Sizes==M[1])]]

  for(i in 1:(ncol(MAT)-1)){
    for(j in (R[i,1]+1):(R[i,2]+1)){
      MAT[j:(j+M[i+1]),i+1] <- MAT[j:(j+M[i+1]),i+1] +MAT[j,i]*probs[[which(Sizes==M[i+1])]]
    }
  }

  return(MAT)
}



#' Discord_Go
#'
#' @param go go
#' @param n n
#' @param fgo fgo
#' @param N N
#' @param prob prob
#' @param Weights weights
#'
#' @return These functions are used by other functions
#' @export
#'
Discord_Go <- function(go,n,fgo,N,prob,Weights=NA){
  if(n==N | is.infinite(go)){
    return(0)
  }

  if(any(is.na(Weights))){
    Weights = dbinom(go:n,n,prob)
    Weights = Weights/sum(Weights)
    warning('NAs in Weights. Weights calculated.')
  }else{
    if(length(Weights)!=n-go+1){
      stop('Weights must be NA or have length n-go+1')
    }
  }

  Prob = Weights *(1-pbinom(fgo-go:n-1,N-n,prob))

  return(1-sum(Prob))
}



#' Discord_Grey
#'
#' @param low low
#' @param up up
#' @param n n
#' @param fgo fgo
#' @param fnogo gnogo
#' @param N N
#' @param prob prob
#' @param Weights weights
#'
#' @return These functions are used by other functions
#' @export
#'
Discord_Grey <- function(low,up,n,fgo,fnogo,N,prob,Weights=NA){
  if(n==N){
    return(0)
  }
  if(is.infinite(low)){
    low <- 0
  }
  if(is.infinite(up)){
    up <- n
  }

  if(any(is.na(Weights))){
    Weights = dbinom(low:up,n,prob)
    Weights = Weights/sum(Weights)
    warning('NAs present in Weights. Weights calculated.')
  }else{
    if(length(Weights)!=up-low+1){
      stop('Weights must be NA or have length up-low+1')
    }
  }
  # browser()

  Prob = Weights *(pbinom(fnogo-low:up,N-n,prob)+1-pbinom(fgo-low:up-1,N-n,prob))

  return(sum(Prob))
}


#' Discord_NoGo
#'
#' @param nogo nogo
#' @param n n
#' @param fnogo fnogo
#' @param N N
#' @param prob prob
#' @param Weights Weights
#'
#' @return These functions are used by other functions
#' @export
#'
Discord_NoGo <- function(nogo,n,fnogo,N,prob,Weights=NA){
  if(n==N | is.infinite(nogo)){
    return(0)
  }

  if(any(is.na(Weights))){
    Weights = dbinom(0:nogo,n,prob)
    Weights = Weights/sum(Weights)
    warning('NAs present in weights. Weights calculated.')
  }else{
    if(length(Weights)!=nogo+1){
      stop('Weights must be NA or have length nogo+1')
    }
  }

  Prob = Weights *pbinom(fnogo-0:nogo,N-n,prob)

  return(1-sum(Prob))
}

