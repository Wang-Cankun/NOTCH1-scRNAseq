#########
# Load useful functions, do not print in the final report
#########
quiet <- function(x) {
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
# point size function from test datasets
x <- c(0,90,124,317,1000,2368,3005,4816,8298,50000,500000,5000000)
y <- c(1,1,0.89,0.33,0.30,0.25,0.235,0.205,0.18,0.1,0.1,0.1)
get_point_size <- approxfun(x, y)

##############################

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return (x)
}

enrichment_dotplot <- function(terms, filename, width=4000, height=2000) {
  terms[, 1] <-
    str_replace_all(terms[, 1], " \\(GO.*", "")
  terms <- terms %>%
    dplyr::filter(Adjusted.P.value < 0.1)
  if (nrow(terms) > 0) {
    new_df <- terms %>%
      mutate(
        Term = as_factor(str_replace_all(Term, " \\(GO.*", "")),
        len = lengths(str_split(Genes, ";")),
        pval = -log10(Adjusted.P.value)
      ) %>%
      rowwise() %>%
      mutate(gene_ratio = eval(parse(text = str_remove_all(Overlap, " ")))) %>%
      dplyr::select(Term, len, pval, gene_ratio, Adjusted.P.value) %>%
      mutate(Term = fct_reorder(Term, Adjusted.P.value))
    
    new_df$Term <-
      factor(new_df$Term, levels = rev(levels(factor(new_df$Term))))
    
    p1 <- ggplot(new_df,
                 aes(x = gene_ratio,
                     y = Term)) +
      geom_point(aes(size = len, color = Adjusted.P.value)) +
      scale_color_gradient(low = "blue",
                           high = "red",
                           trans = 'reverse') +
      theme_bw() +
      ylab("") +
      labs(size = "Overlapping count", color = "Adjusted p-value") +
      theme(
        legend.title = element_text(size = 14, ),
        legend.text = element_text(size = 10, ),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
      ) +
      scale_x_continuous(name = "Overlapping ratio") +
      scale_size(range = c(6, 14))
    
    png(
      filename,
      width = width,
      height = height,
      res = 300
    )
    print(p1)
    dev.off()
  }
}
