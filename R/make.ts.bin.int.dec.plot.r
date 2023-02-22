#' @title Make two-sample binary interim decision plot
#'
#' @param for.plot the results from get.ts.bin.decision
#' @return a ggplot object showing the decision points for the interim analysis.
#' @export
#'
#' @examples
#' my.ts.bin.int.dec <- get.ts.bin.int.dec(
#' a.con = 1, b.con = 1, a.trt = 1, b.trt = 1,
#' Delta.lrv = 0.15, Delta.tv = .30,
#' tau.tv = 0.10, tau.lrv = 0.80, tau.ng = 0.65,
#' go.thresh=.8, ng.thresh=.8,
#' n.trt = 40, n.con = 40,
#' n.int.c = c(10, 20, 30), n.int.t = c(10, 20, 30),
#' runs=500, include_nogo = TRUE)
#' my.ts.bin.int.dec.plot <- make.ts.bin.int.dec.plot(for.plot=my.ts.bin.int.dec)
#' my.ts.bin.int.dec.plot
make.ts.bin.int.dec.plot <- function(for.plot){

        Frm <- for.plot%>%
                mutate(Interim = factor(Interim, levels = 1:length(n.int.c),
                                        labels = paste0('Control = ',n.int.c,', Treatment = ',n.int.t)) )

        levels(Frm$Decision) <- c("Accelerate", "Wait", "Do not Accelerate")
        ggplot(Frm,aes(x=IntermR_C, y = IntermR_T,fill = Decision))+
                geom_tile(aes(color =Decision))+
                facet_wrap(~Interim,scales = 'free')+
                scale_y_continuous(expand = c(0,0), breaks=pretty(0:max(Frm$IntermR_C)), minor_breaks = NULL)+
                scale_x_continuous(expand = c(0,0), breaks=pretty(0:max(Frm$IntermR_T)), minor_breaks = NULL)+
                # scale_y_continuous(expand = c(0,0),breaks = function(limits){(limits[1]+.5):(limits[2]-.5)}, minor_breaks = NULL)+
                # scale_x_continuous(expand = c(0,0),breaks = function(limits){(limits[1]+.5):(limits[2]-.5)}, minor_breaks = NULL)+
                scale_fill_manual(values = c(alpha("green", .5),alpha("grey", .5),alpha("red", .5)))+
                scale_color_manual(values = c('darkgreen','grey70','pink'))+
                xlab('Successes in Control')+
                ylab('Successes in Treatment')+
                theme_bw()+
                theme(legend.position = 'bottom')
}
