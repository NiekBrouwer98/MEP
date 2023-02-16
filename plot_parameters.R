
point <- ggplot(success_models[[1]][[3]]) +
  geom_point(aes(x=a, y=b, color=phenotype_combo), alpha=0.7, size=0.5) +
  # geom_density_2d() +
  theme_bw() + scale_y_log10() +
  xlab('Shape') + ylab('Scale') + scale_color_brewer(palette='Set2') +
  theme(legend.position='right') + ggtitle(c) +ylim(0,300) + xlim(0,5)

print(point)
