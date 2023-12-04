#!/bin/bash

#'
#' Calculate repeat divergences.
#'


export REPEATS_IN_DIR=/staging/lnell/annotations
export OUT_DIR=repeats_diverg


export REPEATS_TARS=$(find ${REPEATS_IN_DIR} -name *_repeats.tar.gz)

mkdir ${OUT_DIR}
cd ${OUT_DIR}

for rep_tar in ${REPEATS_TARS}; do
    rt_dir=$(basename ${rep_tar%.tar.gz})
    tar -xf ${rep_tar} -C ./
    check_exit_status "null" $?
    mv ${rt_dir}/masker/*.align ./${rt_dir/_repeats/}.align
    check_exit_status "null" $?
    rm -r ${rt_dir}
    check_exit_status "null" $?
    # (interactive shell)
    if [ "$out_loc" != "/dev/null" ]; then echo $rt_dir; fi
    unset -v rt_dir
done

. /app/.bashrc
conda activate repeat-env
export PERL5LIB=/opt/conda/envs/repeat-env/share/RepeatMasker

# Based on my naming scheme, this should return the species names:
export SPP_NAMES=$(find . -type f -name "*.align" | sed 's/\.\///g;s/\.align//g')


for spp in ${SPP_NAMES}; do
    calcDivergenceFromAlign.pl -s ${spp}.divsum ${spp}.align
done



# Pretend python dictionary of genome sizes (not including Ns):
export gsizes=("Aaegyp:1278709169" \
               "Asteph:243467213" \
               "Bantar:88894437" \
               "Cmarin:84365937" \
               "Cquinq:547064555" \
               "Cripar:191805449" \
               "Csonor:155941109" \
               "Ctenta:179550640" \
               "Mdomes:691720152" \
               "Pakamu:85833737" \
               "Ppemba:122659044" \
               "Pstein:143569098" \
               "Pvande:118352072" \
               "Tgraci:91827299")

for gs in ${gsizes[@]}; do
    spp=${gs%%:*}
    gsize=${gs#*:}
    createRepeatLandscape.pl -div ${spp}.divsum -g $gsize \
        > ${spp}.html
    check_exit_status "$spp landscape" $?
done



rm *.align

mkdir divsum && mv *.divsum ./divsum/

mkdir html && mv *.html ./html/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz "${REPEATS_IN_DIR}/"

rm -r ${OUT_DIR}


