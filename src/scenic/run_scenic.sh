#!/bin/bash
set -e -o pipefail

seu_path=`readlink -f $1`

script_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
output_dir=$( cd "$(dirname $(dirname "$seu_path"))" ; pwd -P )

echo $parent_path
echo $seu_path
echo $output_dir

dataset=`basename $seu_path .rds`

input_loom_path=`realpath $output_dir/scenic/$dataset.loom`

output_loom=`echo $dataset-final.loom`
output_loom_path=`realpath $output_dir/scenic/$output_loom`


Rscript "$script_path"/convert_seu_to_loom_scenic.R $seu_path

echo "$input_loom_path"
echo "$output_loom_path"

# nextflow run aertslab/SCENICprotocol \
#         -profile docker \
#         --loom_input "$input_loom_path" \
#         --loom_output "$output_loom" \
#         --TFs ~/Homo_sapiens/scenic/hs_hgnc_curated_tfs.txt \
#         --motifs ~/Homo_sapiens/scenic/motifs-v9-nr.hgnc.tbl \
#         --db ~/Homo_sapiens/scenic/hg38_refseq-10kb_up_and_down_tss.mc9nr.feather \
#         --thr_min_genes 1

# for testing file paths and script
nextflow run aertslab/SCENICprotocol \
    -profile docker \
    --loom_input src/scenic_src/example/expr_mat_tiny.loom \
    --loom_output "$output_loom" \
    --TFs src/scenic_src/example/test_TFs_tiny.txt \
    --motifs src/scenic_src/example/motifs.tbl \
    --db src/scenic_src/example/*feather \
    --thr_min_genes 1

tmux_session=`tmux display-message -p "#S"`
echo $tmux_session
if [ $? == 1 ]
then
    echo -e "Subject:failure\n$tmux_sesssion" | sendmail -t kevin.stachelek@gmail.com
else
    echo -e "Subject:success!\n$tmux_session" | sendmail -t kevin.stachelek@gmail.com
fi


