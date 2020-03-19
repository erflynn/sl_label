#!/bin/bash
#SBATCH --job-name=convert_chunk_%A_%a
#SBATCH --output=logs/convert_chunk_%A_%a.out
#SBATCH --error=logs/convert_chunk_%A_%a.err
#SBATCH --time=6:00:00 
#SBATCH --mem=40000
#SBATCH --partition=rbaltman
#
# Grab a chunk of the input .tsv data and put it into a .gct file with the 
# proper header setup.


infile=$1
prefix=$2
inc=$3
i=$SLURM_ARRAY_TASK_ID


# get the chunk start/end
ncol=$(awk 'NR==3{print NF}' "$infile")
nrow=$(awk 'END{print NR}' "$infile")
start=$((inc*i+1))
end=$((inc*(i+1)+1))
if [ "$end" -gt "$ncol" ]; then
    end=$ncol;
fi

# read in a chunk
echo $start;
echo $end;
cut -f1,$start-$end "$infile" > "${infile}.$i";

echo "CUT";

my_col=$((end-start-1))
my_row=$((nrow-1))
head -1 ${infile}.$i | awk '{printf "GID%s\n", $0}' > ${prefix}_tmp_${i}.txt;
cat <(echo -e "#1.2\n${my_row}\t${my_col}") ${prefix}_tmp_${i}.txt > ${prefix}_cols_${i}.txt;
rm ${prefix}_tmp_${i}.txt

echo "header created"
# add header for gct

tail -n +2 ${infile}.$i > ${prefix}_${i}v2.gct;

cat ${prefix}_cols_${i}.txt ${prefix}_${i}v2.gct > ${prefix}_${i}.gct;
