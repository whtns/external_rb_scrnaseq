# Four-panel figure supporting the numbat round-selection policy.
#   A  segments called per round (paired, per sample)
#   B  round-1 segments that survive to round 2 vs those that do not
#   C  consecutive-round segment-set agreement - the convergence test
#   D  the canonical RB events a fixed i=2 misses, by arm coverage and LLR
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

BY  <- fread("results/round_convergence_by_round.csv")
P   <- fread("results/round1_segment_persistence.csv")
S   <- fread("results/round_convergence_stability.csv")
M   <- fread("results/round_policy_missed_events.csv")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

## A -- segments per round
a <- ggplot(BY, aes(factor(round), n_segs)) +
  geom_line(aes(group = sample_id), colour = "grey70", linewidth = .3) +
  geom_point(colour = "grey45", size = .7) +
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               colour = "firebrick", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", colour = "firebrick", size = 2.4) +
  labs(title = "A  segments called per round",
       subtitle = "non-neutral consensus segments; red = mean",
       x = "numbat round", y = "segments") + th

## B -- fate of round-1 segments
P[, fate := ifelse(persists, "kept at r2", "dropped at r2")]
b1 <- ggplot(P, aes(fate, size_mb, fill = fate)) +
  geom_boxplot(outlier.size = .4, width = .55, show.legend = FALSE) +
  scale_y_log10() +
  scale_fill_manual(values = c("dropped at r2" = "#d98a8a", "kept at r2" = "#8aa9d9")) +
  labs(x = NULL, y = "segment size (Mb, log)") + th
b2 <- ggplot(P, aes(fate, LLR, fill = fate)) +
  geom_boxplot(outlier.size = .4, width = .55, show.legend = FALSE) +
  scale_y_log10() +
  scale_fill_manual(values = c("dropped at r2" = "#d98a8a", "kept at r2" = "#8aa9d9")) +
  labs(x = NULL, y = "LLR (log)") + th
b <- (b1 | b2) + plot_annotation(title = "B  fate of round-1 segments")

## C -- convergence test
c1 <- ggplot(S, aes(pair, jaccard)) +
  geom_line(aes(group = sample_id), colour = "grey80", linewidth = .3) +
  geom_boxplot(width = .45, fill = "#9ecae1", outlier.size = .5, alpha = .8) +
  stat_summary(fun = mean, geom = "point", colour = "firebrick", size = 2.4) +
  ylim(0, 1) +
  labs(title = "C  agreement between consecutive rounds",
       subtitle = "no plateau: still rising at r3->r4 (p = 0.039)",
       x = NULL, y = "Jaccard (>=50% reciprocal overlap)") + th

## D -- what a fixed i=2 misses
M[, lab := paste0(sample_id, "  ", event)]
setorder(M, status, arm_frac_rk)
M[, lab := factor(lab, levels = lab)]
d <- ggplot(M, aes(arm_frac_rk, lab, colour = status)) +
  geom_segment(aes(x = arm_frac_r2, xend = arm_frac_rk, yend = lab),
               colour = "grey75", linewidth = .4) +
  geom_point(aes(x = arm_frac_r2), colour = "grey45", size = 1.1) +
  geom_point(size = 1.9) +
  geom_vline(xintercept = 0.15, linetype = "dashed", colour = "firebrick") +
  scale_colour_manual(values = c("absent at r2" = "#c0392b",
                                 "sub-threshold at r2" = "#e08214")) +
  labs(title = "D  canonical RB events a fixed i = 2 would miss",
       subtitle = "grey = coverage at r2, colour = at r1; dashed = 0.15 arm floor",
       x = "fraction of target arm covered", y = NULL, colour = NULL) +
  th + theme(legend.position = "bottom", axis.text.y = element_text(size = 5.5))

pdf("results/round_policy_report.pdf", width = 11, height = 13)
print((a | c1) / wrap_elements(b) / d + plot_layout(heights = c(1, 1, 2.1)))
dev.off()
cat("wrote results/round_policy_report.pdf\n")
