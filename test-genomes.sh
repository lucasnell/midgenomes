

mkdir working
cd working


# cp /staging/lnell/Cson.tar.gz ./
# tar -xzf Cson.tar.gz
# rm Cson.tar.gz
#
# # Rename fastas to be more descriptive for downstream naming:
# mv Culicoides_sonorensis.Cson1.dna_sm.toplevel.fa ensemble.fasta
# mv GCA_900258525.3_C_sonorensis_v2_redundans_genomic.fa genbank_v2.fasta
# mv GCA_900002565.1_ASM90000256v1_genomic.fa genbank_v1.fasta

cp /staging/lnell/Cripa.tar.gz ./
tar -xzf Cripa.tar.gz
mv Cripa/Cripa*.fasta ./
rm -r Cripa.tar.gz Cripa


export n_threads=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')
# export n_threads=4

# conda activate busco-env

busco --download diptera_odb10


for fasta in *.fasta; do
    # summ-scaffs.py ${fasta} | tee summ_${fasta/.fasta/}.out
    # check_exit_status "summ-scaffs.py ${fasta}" $?
    busco -m genome -l ./busco_downloads/lineages/diptera_odb10 -i ${fasta} \
        -o busco_${fasta/.fasta/} --cpu ${n_threads} --offline | \
        tee busco_${fasta/.fasta/}.out
    check_exit_status "busco ${fasta}" $?
done

