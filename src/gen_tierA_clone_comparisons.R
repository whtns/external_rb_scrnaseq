# Generate config entries so the tier-A (tumor-vs-tumor) RB SCNA contrasts get
# two-clone collages and diffex, alongside the existing curated ones.
#
# Emits, without overwriting anything in place:
#   results/diploid_audit/large_clone_comparisons_MERGED.yaml
#   results/diploid_audit/scna_collage_samples_PROPOSED.R
# Existing YAML entries always win: this only ADDS comparison keys a sample does
# not already have, so the curated pairs and their current collages never move.

suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr); library(tidyr)})
setwd("/project2/cobrinik_1090/external_rb_scrnaseq_proj")
out_dir <- "results/diploid_audit"

RB  <- c("1q+","2p+","6p+","16q-")
key <- c("1q+"="1q","2p+"="2p","6p+"="6p","16q-"="16q")
centromeres_bp <- c("1"=123500000L,"2"=93100000L,"3"=92200000L,"4"=50700000L,
  "5"=48300000L,"6"=59200000L,"7"=59800000L,"8"=45000000L,"9"=44400000L,
  "10"=40600000L,"11"=52700000L,"12"=35600000L,"13"=17000000L,"14"=17100000L,
  "15"=18400000L,"16"=37300000L,"17"=24900000L,"18"=18200000L,"19"=25900000L,
  "20"=28200000L,"21"=11900000L,"22"=14000000L)
lab_of <- function(chrom,start,end,st){
  chr <- as.character(chrom); mid <- (start+end)/2
  arm <- ifelse(!chr %in% names(centromeres_bp),"?",ifelse(mid<centromeres_bp[chr],"p","q"))
  sfx <- case_when(st %in% c("amp","bamp")~"+", st %in% c("del","bdel")~"-",
                   st %in% c("loh","cnloh")~"cnloh", TRUE~st)
  paste0(chr,arm,sfx)
}

sel <- read_csv("results/numbat_selected_round.csv", show_col_types=FALSE) |>
  select(sample_id, selected_round)
elig <- read_csv(file.path(out_dir,"rb_scna_comparison_eligibility.csv"), show_col_types=FALSE) |>
  filter(tier == "A")

# recover the seg_cons tokens each contrast gains for its RB SCNA
tok_for <- function(sid, rnd, cw, cwo, rb) {
  cp <- read.delim(sprintf("output/numbat_sridhar/%s/clone_post_%s.tsv", sid, rnd),
                   stringsAsFactors=FALSE, colClasses="character")
  sc <- read.delim(sprintf("output/numbat_sridhar/%s/segs_consensus_%s.tsv", sid, rnd),
                   stringsAsFactors=FALSE)
  cp$GT_opt[is.na(cp$GT_opt)] <- ""
  cp$clone_opt <- as.integer(cp$clone_opt)
  g <- function(c) { v <- unique(cp$GT_opt[cp$clone_opt == c]); if (!length(v)||v[1]=="") character(0) else str_split_1(v[1],",") }
  gained <- setdiff(g(cw), g(cwo))
  sn <- sc[!is.na(sc$cnv_state) & sc$cnv_state != "neu", , drop=FALSE]
  m  <- setNames(lab_of(sn$CHROM, sn$seg_start, sn$seg_end, sn$cnv_state), sn$seg_cons)
  gained[!is.na(m[gained]) & m[gained] == rb]
}

elig$tokens <- mapply(function(s,w,wo,r){
  rnd <- sel$selected_round[sel$sample_id==s][1]
  paste(tok_for(s,rnd,w,wo,r), collapse=",")
}, elig$sample_id, elig$clone_w, elig$clone_wo, elig$rb_gained)

elig$comp_key <- sprintf("%s_v_%s_%s", elig$clone_w, elig$clone_wo, elig$rb_gained)

y <- yaml::read_yaml("config/large_clone_comparisons.yaml")
added <- 0L; skipped <- 0L; notok <- 0L
for (i in seq_len(nrow(elig))) {
  s <- elig$sample_id[i]; k <- elig$comp_key[i]
  toks <- if (nzchar(elig$tokens[i])) str_split_1(elig$tokens[i], ",") else character(0)
  if (!length(toks)) { notok <- notok + 1L
    message("no tokens resolved for ", s, " ", k, " -- skipped"); next }
  if (is.null(y[[s]])) y[[s]] <- list()
  if (!is.null(y[[s]][[k]])) { skipped <- skipped + 1L; next }   # curated entry wins
  y[[s]][[k]] <- as.list(toks)
  added <- added + 1L
}
writeLines(yaml::as.yaml(y), file.path(out_dir,"large_clone_comparisons_MERGED.yaml"))

cat("tier A contrasts:", nrow(elig), "\n added:", added,
    " already present:", skipped, " unresolved tokens:", notok, "\n\n")
print(as.data.frame(elig |> select(sample_id, scna, comp_key, tokens, n_cells_wo, n_cells_w)), row.names=FALSE)

# proposed collage sample sets = existing curated set U tier A
cur <- list(
  "1q"=c("SRX10264523","SRX10264526","SRX11133594","SRX11133593","SRX11133592","SRX10831287"),
  "2p"=c("SRX10031193","SRX10264517","SRX10264518","SRX10264519","SRX10264520","SRX10264523","SRX14116946","SRX14116947","SRX22868102"),
  "6p"=c("SRX10031193","SRX10831281","SRX10831282","SRX11133588","SRX14116944","SRX14116946","SRX14116947","SRX22868105"),
  "16q"=c("SRX11133594","SRX11133593","SRX11133592"))
prop <- lapply(names(cur), function(k)
  sort(union(cur[[k]], elig$sample_id[elig$scna == k])))
names(prop) <- names(cur)
lines <- c("    tar_target(scna_collage_samples,", "      list(")
for (k in names(prop)) {
  ids <- paste0('"', prop[[k]], '"', collapse = ", ")
  lines <- c(lines, sprintf('        "%s" = c(%s)%s', k, ids,
                            if (k == tail(names(prop),1)) "" else ","))
}
lines <- c(lines, "      )", "    ),")
writeLines(lines, file.path(out_dir,"scna_collage_samples_PROPOSED.R"))
cat("\n=== proposed collage set sizes ===\n")
for (k in names(prop)) cat(sprintf("  %-4s %2d -> %2d\n", k, length(cur[[k]]), length(prop[[k]])))
