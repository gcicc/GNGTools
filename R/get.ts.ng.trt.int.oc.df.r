#' @title Get two-sample normal-gamma interim OC curve data.frame
#'
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size parameter for treatment group
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size parameter for control group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param mu.c assumed mean for control group
#' @param s.t treatment standard deviation
#' @param s.c control standard deviation
#' @param npointsLookup number of points for lookup table
#' @param npoints number of points to run simulations
#' @param n.MC.lookup number of trials used for lookup table
#' @param n.MC number of trials run at each point
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param n.int.t interim sample sizes for treatment arm
#' @param n.int.c interim sample sizes for control arm
#' @param final.n.t final sample size: treatment arm
#' @param final.n.c final sample sizE: control arm
#' @param go.thresh interim predictive probability threshold for go
#' @param ng.thresh interim predictive probability threshold for no-go
#' @param go.parallel logical for parallel processing
#' @param cl cl
#' @param include_nogo logical
#' @param seed random seed
#'
#' @return A data.frame is returned.
#' @export
#'
#' @examples \donttest{
#' my.ts.ng.trt.int.oc.df <- get.ts.ng.trt.int.oc.df(npointsLookup = 2, npoints=3, n.MC.lookup=5,
#' n.MC=5, go.parallel=FALSE)
#' }
get.ts.ng.trt.int.oc.df <- function(
                mu.0.t = 0, n.0.t = .0001, alpha.0.t = 0.25,beta.0.t = 1,
                mu.0.c = 0, n.0.c = 10, alpha.0.c = 2.5, beta.0.c = 10,
                Delta.lrv = 1.5, Delta.tv = 3.0,
                mu.c=0.25, s.t = 1.5, s.c = 1.5,
                npointsLookup = 20, npoints=20, n.MC.lookup=500, n.MC=500,
                tau.tv = 0.1, tau.lrv = 0.8,tau.ng = 0.65,
                n.int.t = c(27, 40), n.int.c = c(27, 40),
                final.n.t = 55, final.n.c = 55,
                go.thresh = 0.6,ng.thresh = 0.6,
                go.parallel = TRUE, cl=cl,
                seed=1234, include_nogo=TRUE)
{
        set.seed(seed)
        message("0. Create Look-up table\n") #----
        # 0. Create a lookup table of studyend critera over a grid of study.end control means
        #

                  #tictoc::tic() #  47.34 sec elapsed with 20 = npoints.lookup
        # This creates a lookup grid hold requirements for study-end xbar.t for a grid of study.end control means.
        # In simulations we appeal to the lookup grid, find the closest xbar.c in the lookup table and pick off the required studyend go and no go
        look.up.table <- make.ts.ng.studyend.criteria(
                mu.0.c = mu.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c, n.0.c = n.0.c,
                mu.0.t = mu.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                xbar.c.LB = mu.c - 5*s.c/sqrt(n.int.c[1]),
                xbar.c.UB = mu.c + 5*s.c/sqrt(n.int.c[1]),
                npoints=npointsLookup, s.c = s.c, n.c = final.n.c,
                xbar.t = mu.0.t, s.t = s.t, n.t = final.n.t,
                Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                tau.ng = tau.ng, tau.lrv = tau.lrv, tau.tv = tau.tv,n.MC = n.MC.lookup,
                go.parallel=go.parallel, cl=cl)


                  message("1. Create grid of control and treatment effects\n") # ----
                  # Grid size is controlled by npoints and grows quadratically!
                  # Go from control mean to control mean + 1.5*Delta.tv
                  # Repeat, this time starting from mu.c-s.c
                  data.in <- bind_rows(
                          tibble(mu.0.c = mu.0.c, alpha.c = alpha.0.c, beta.c = beta.0.c, n.0.c = n.0.c,
                                 mu.0.t = mu.0.t, alpha.t = alpha.0.t, beta.t = beta.0.t, n.0.t = n.0.t,
                                 s.t=s.t, s.c=s.c, final.n.t=final.n.t, final.n.c=final.n.c,
                                 expand.grid(mu.c=mu.c, mu.t=seq(mu.c, mu.c + Delta.tv*1.5, length.out=npoints))),
                          tibble(
                                  mu.0.c = mu.0.c, alpha.c = alpha.0.c, beta.c = beta.0.c, n.0.c = n.0.c,
                                  mu.0.t = mu.0.t, alpha.t = alpha.0.t, beta.t = beta.0.t, n.0.t = n.0.t,
                                  s.t=s.t, s.c=s.c,final.n.t=final.n.t, final.n.c=final.n.c,
                                  expand.grid(mu.c=mu.c-s.c, mu.t=seq(mu.c-s.c, mu.c-s.c + Delta.tv*1.5, length.out=npoints))),
                          tibble(mu.0.c = mu.0.c, alpha.c = alpha.0.c, beta.c = beta.0.c, n.0.c = n.0.c,
                                 mu.0.t = mu.0.t, alpha.t = alpha.0.t, beta.t = beta.0.t, n.0.t = n.0.t,
                                 s.t=s.t, s.c=s.c, final.n.t=final.n.t, final.n.c=final.n.c,
                                 expand.grid(mu.c=mu.c+s.c, mu.t=seq(mu.c+s.c, mu.c+s.c + Delta.tv*1.5, length.out=npoints))))

        message("2. Build trial data\n") # -----
        # 2. Build trial data
        # For each row of data.in (our simulation grid), we run n.MC trials
        message("2a. Define function\n") # ----
        # -------function----
        get.interim.results <- function(IA.c = c(20, 30, 40),
                                        IA.t = c(20, 30, 40),
                                        mu.c=0, mu.t=1,
                                        sigma.c = 2, sigma.t=2){
                # CHECK: length(IA.c) == length(IA.t)
                interim.results <- data.frame(t(apply(matrix(1:length(IA.c)), MARGIN=1, function(x){
                        # Get full.data
                        temp.c <- rnorm(n = IA.c[length(IA.c)], mean=mu.c, sd=sigma.c)
                        temp.t <- rnorm(n = IA.t[length(IA.t)], mean=mu.t, sd=sigma.t)
                        matrix(data = c(
                                n.c=IA.c[x],
                                mean.c=mean(temp.c[1:IA.c[x]]),
                                sd.c=sd(temp.c[1:IA.c[x]]),
                                n.t=IA.t[x],
                                mean.t=mean(temp.t[1:IA.t[x]]),
                                sd.t=sd(temp.t[1:IA.t[x]])),nrow=1, ncol=6)
                        })))
                colnames(interim.results) <- c("int.n.c", "xbar.c", "sd.c", "int.n.t", "xbar.t", "sd.t")
                interim.results <- interim.results
                return(interim.results)
        }

        message("2b. Create trials\n") # ----
        # npoints = 5, n.MC = 500: 5.61 sec elapsed sec elapsed with 3 interims + final dim(my.trials) = 200000 x 7
        # npoints = 10, n.MC = 500: 21.92 sec elapsed sec elapsed with 3 interims + final dim(my.trials) = 200000 x 7
        # npoints = 10, n.MC = 1000: 52.77 sec elapsed with 3 interims + final dim(my.trials) = 400000 x 7

        if(go.parallel==F){
                #tictoc::tic() # 3.69 sec elapsed
                my.trials <- bind_rows(
                        apply(X = matrix(1:nrow(data.in)), MARGIN=1, FUN = function(y){
                                data.frame(do.call(rbind,
                                                   (lapply(X=1:n.MC, FUN = function(x) {
                                                           get.interim.results(
                                                                   IA.c = c(n.int.c, final.n.c),
                                                                   IA.t = c(n.int.t, final.n.t),
                                                                   mu.c=data.in$mu.c[y], mu.t=data.in$mu.t[y],
                                                                   sigma.c = data.in$s.t[y], sigma.t=data.in$s.t[y])})))) %>%
                                        mutate(mu.c = data.in$mu.c[y], sigma.c=data.in$s.t[y], mu.t=data.in$mu.t[y], sigma.t=data.in$s.t[y])

                                      })) %>%
                        mutate(trial=rep(1:(nrow(.)/length(c(n.int.c, final.n.c))), each=length(c(n.int.c, final.n.c))))
                #tictoc::toc()
                } else {
                        #tictoc::tic() # 1.02 sec elapsed
                        clusterEvalQ(cl = cl, expr = {requireNamespace("dplyr"); requireNamespace("tidyr")})
                        clusterExport(cl = cl, varlist = c("get.ts.ng.studyend.GNG", "get.ts.ng.mc.df","get.ts.ng.mc",
                                                           "get.ng.post","rt_ls","get.interim.results",
                                                           "mu.0.c" , "alpha.0.c", "beta.0.c", "n.0.c",
                                                           "mu.0.t", "alpha.0.t", "beta.0.t", "n.0.t",
                                                           "s.c", "s.t", "Delta.lrv" , "Delta.tv",
                                                           "tau.ng", "tau.lrv", "tau.tv",
                                                           "n.MC", "n.int.c", "n.int.t", "final.n.c",
                                                           "final.n.t", "data.in"), envir = environment())
                        my.trials <- bind_rows(
                                parApply(cl=cl, X = matrix(1:nrow(data.in)), MARGIN=1, FUN = function(y){
                                        data.frame(do.call(rbind,
                                                           (lapply(X=1:n.MC, FUN = function(x) {
                                                                   get.interim.results(
                                                                           IA.c = c(n.int.c, final.n.c),
                                                                           IA.t = c(n.int.t, final.n.t),
                                                                           mu.c=data.in$mu.c[y], mu.t=data.in$mu.t[y],
                                                                           sigma.c = data.in$s.t[y], sigma.t=data.in$s.t[y])})))) %>%
                                                mutate(mu.c = data.in$mu.c[y], sigma.c=data.in$s.t[y], mu.t=data.in$mu.t[y], sigma.t=data.in$s.t[y])
                                        })) %>%
                        mutate(trial=rep(1:(nrow(.)/length(c(n.int.c, final.n.c))), each=length(c(n.int.c, final.n.c))))
                        #tictoc::toc()
                        }

  # Get study end results
  # Accomplish fast merge of look-up table with trial.results

  message("3. Get study end results and employ lookup table\n")
  #tictoc::tic() # 3.33 sec elapsed
  study.end.results <- data.table(my.trials %>% group_by(trial) %>% slice(n())) %>% dplyr::select(trial, xbar.c, xbar.t)

  lookup.go <- data.table(look.up.table %>% dplyr::filter(result=="Go") %>% dplyr::select(xbar.c, xbar.t))
  lookup.nogo <- data.table(look.up.table %>% dplyr::filter(result=="No-Go") %>%  dplyr::select(xbar.c, xbar.t))
  lookup.go[,merge:=xbar.c]
  lookup.nogo[,merge:=xbar.c]
  study.end.results[,merge:=xbar.c]
  setkeyv(study.end.results,c('merge'))
  setkeyv(lookup.go,c('merge'))
  setkeyv(lookup.nogo,c('merge'))
  Merge_GO=lookup.go[study.end.results ,roll='nearest']
  Merge_NOGO=lookup.nogo[study.end.results,roll='nearest']
  # find the closesest match here... and pickoff the corresponding value of xbar.t

  my.trials <- my.trials %>%
    left_join(Merge_NOGO %>% dplyr::select(trial, xbar.t) %>% dplyr::rename(xbar.NG=xbar.t)) %>%
    left_join(Merge_GO %>% dplyr::select(trial, xbar.t) %>% dplyr::rename(xbar.GO=xbar.t)) %>%
          arrange(trial, int.n.c)
  #tictoc::toc()


  message("4. Process sims\n")

  # Check on study end
  # my.trials %>% dplyr::filter(int.n.c==55) %>% ggplot(aes(x=xbar.c, y=xbar.t, color=xbar.t > xbar.GO))+ geom_point()
  # NOTE mu.c plays a role...
  # process.sim %>% group_by(mu.c) %>% summarize(ng = mean(xbar_ng), go=mean(xbar_go)) %>% mutate(treat.ng= ng-mu.c, treat.go = go-mu.c)

  process.sim <- my.trials %>%
    dplyr::select(trial, mu.c, sigma.c, mu.t, sigma.t, xbar.c, sd.c, xbar.t, sd.t,  xbar.NG, xbar.GO) %>%
    bind_cols(
      return.ts.ng.int.req.df(mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
                         mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c,
                         xbar.t=my.trials$xbar.t,
                         s.t = my.trials$sd.t,
                         n.t =   my.trials$int.n.t,
                         xbar.c=my.trials$xbar.c,
                         s.c = my.trials$sd.c,
                         n.c =  my.trials$int.n.t,
                         Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                         tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                         xbar_ng = my.trials$xbar.NG,
                         xbar_go = my.trials$xbar.GO,
                         go.thresh=go.thresh, ng.thresh=ng.thresh,
                         n.MC = n.MC) %>%
        dplyr::select(mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, n.c, n.t, Go, NoGo, decision)) %>%
    dplyr::select( trial, # trial
                   mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, # Prior specs
                   mu.c, sigma.c, mu.t, sigma.t, # data generation specs
                   xbar.c, sd.c, n.c, xbar.t, sd.t,  n.t, xbar.NG,  xbar.GO, Go, NoGo, decision)
  #tictoc::toc()
  message("5. Prep summary\n")

  #tictoc::tic() # 7.84 sec elapsed
  # This reports any analysis with priority to first Go/No-Go among interims and final
  any.analysis <- bind_rows(
    # This batch takes the first non-consider
    process.sim %>%
            group_by(trial) %>%
            mutate(decision = ifelse(n.t==max(n.t),
                                     case_when(
                                             xbar.t < xbar.NG ~ "No-Go",
                                             xbar.t > xbar.GO ~ "Go",
                                             TRUE ~ "Consider"),
                                     decision)) %>%
            group_by(trial, mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
            mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
            mutate(dnum = as.numeric(decision)) %>%
            dplyr::filter(dnum != 3)%>% slice(1),
    # This batch holds the 'all consider' cases
    process.sim %>%
            group_by(trial) %>%
            dplyr::mutate(decision = ifelse(n.t==max(n.t),
                                     case_when(
                                             xbar.t < xbar.NG ~ "No-Go",
                                             xbar.t > xbar.GO ~ "Go",
                                             TRUE ~ "Consider"),
                                     decision)) %>%
            group_by(trial, mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
            dplyr::mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
            dplyr::mutate(dnum = as.numeric(decision)) %>%
            summarize(Consider = all(dnum==3)) %>%
            dplyr::filter(Consider == TRUE) %>%
            dplyr::rename(decision = Consider) %>%
            mutate(decision="Consider"))  %>%
          group_by( mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
          dplyr::count(decision=decision, .drop=FALSE) %>% mutate(analysis="Any Analysis")

      # Broken down by trial with adjustments for previous interims.
  any.interim <-  bind_rows(
          process.sim %>% group_by(trial) %>% slice(-n()) %>%
                  group_by(trial, mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
                  mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
                  mutate(dnum = as.numeric(decision)) %>%
                  dplyr::filter(dnum != 3) %>%
                  dplyr::select(trial, mu.c, sigma.c,  mu.t, sigma.t, decision) %>% slice(1),
          process.sim %>% group_by(trial) %>% slice(-n()) %>%
                  group_by(trial, mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
                  mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
                  mutate(dnum = as.numeric(decision)) %>%
                  dplyr::filter(all(dnum==3))  %>%
                  summarize(Consider = all(dnum==3)) %>%
                  dplyr::filter(Consider == TRUE) %>%
                  dplyr::rename(decision = Consider) %>%
                  dplyr::mutate(decision="Consider")) %>%
          group_by(mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
          dplyr::count(decision=decision, .drop=FALSE) %>%
          dplyr::mutate(analysis="Any Interim")


  # Broken down by interim.n.t, with no adjustments
  each.interim <- process.sim %>% group_by(mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t, n.c) %>%
          dplyr::mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
          dplyr::count(decision=decision, .drop=FALSE) %>%
          dplyr::mutate(analysis=as.character(n.c)) %>%
          ungroup() %>%
          dplyr::select(-n.c)

  # Here's study-end results per simulation
  study.end <-     process.sim %>%
          group_by(trial, mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
          slice(n())  %>%
          dplyr::mutate(decision = case_when(
                  xbar.t < xbar.NG ~ "No-Go",
                  xbar.t > xbar.GO ~ "Go",
                  TRUE ~ "Consider")) %>%
          group_by(mu.0.c, n.0.c, alpha.0.c, beta.0.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, mu.c, sigma.c, mu.t, sigma.t) %>%
          dplyr::count(decision=decision, .drop = FALSE)%>%
          dplyr::mutate(analysis="Study-end")

  for.return <- bind_rows(each.interim, any.interim, study.end, any.analysis) %>%
    group_by(mu.c, sigma.c,  mu.t, sigma.t, analysis) %>%
    mutate(relFreq=n/sum(n),
           Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
           tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

  message("6. return results \n")

  return(for.return)
}
