#' mixture coefficient of measurement distribution
#'
#' @param R states under experiment
#' @param R_0 states in WildType
#' @param flag flag for printing plots
#'
#' @return level of activity
#' @export
gammas = function(R,
                  R_0,
                  flag) {

  if(class(R_0)=="numeric"){
    #make R_0 matrix in case vector being passed
    R_0 <- matrix(R_0,nrow = nrow(R),ncol=ncol(R),byrow = TRUE)
    gamma_val = tanh(abs(R - R_0))

  }else{
    gamma_val = tanh(abs(R - R_0))

  }

  #plot gamma_vals
  if(flag){

    df_states=reshape2:::melt.matrix(gamma_val)

    colnames(df_states) = c("Time","Hidden_node","Activity_Value")
    Hidden_nodes=as.factor(df_states$Hidden_node)
    Times=as.factor(df_states$Time)
    qual_col_pals = RColorBrewer::brewer.pal.info[
      c("Set1","Dark2"),]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal,
                               qual_col_pals$maxcolors,
                               rownames(qual_col_pals)))


    p=ggplot2::ggplot(df_states,
                      ggplot2::aes(x = Times,
                                   y = Activity_Value,
                                   fill=Hidden_nodes,
                                   shape=Hidden_nodes,
                                   color=Hidden_nodes)) +
      ggplot2::geom_point(alpha = 1,
                          size = 2) +
      #theme_minimal() +
      ggplot2::theme(legend.position = "top",
            #aspect.ratio=1,
            axis.text=ggplot2::element_text(size=6),
            axis.title=ggplot2::element_text(size=6),
            plot.title = ggplot2::element_text(size=10),
            legend.title = ggplot2::element_text(size = 6),
            legend.text = ggplot2::element_text(size = 6))+
      ggplot2::ggtitle(paste0("Experiment_",
              which(apply(R, 2, function(x) all(x==0)))))+
      #ggplot2::geom_jitter()+
      ggplot2::scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,
                                             9,10,11,12,13,14,15))+
      ggplot2::scale_color_manual(values=col_vector)+
      ggplot2::scale_fill_manual(values=col_vector)

    print(p)

  }

  return(gamma_val)
}
