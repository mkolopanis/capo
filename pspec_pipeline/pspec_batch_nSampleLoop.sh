#!/bin/bash

### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>

PATH2CAPO='/home/saulkohn/ReposForCanopy/capo/'

#PSA64
if false
then
CALFILE='psa6240_v003'
RA='.1_8.6'
SEP='0,1'
DATA=$1
EVEN_FILES=${DATA}'/even/sep'${SEP}'/*.uvGA'
ODD_FILES=${DATA}'/odd/sep'${SEP}'/*.uvGA'
DIRNAME=$2
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
CHAN='95_115'
#NBOOT=20
POL='I'
weight='I'
WINDOW='none'
FRF=''
#NOFRFPATH='--nofrfpath pspec128_uvGA/inject_sep'${SEP}'_0.01/pspec_pk_k3pk.npz' # path to one pspec_2d_to_1d.py output for NONFRF case
NOFRFPATH=''
LMODE='' #'--lmode=12'
CHANGEC='' #'--changeC'

else
### My Options ###
#CALFILE='psa6622_v003'
#RA='4_10'
#SEP='0,2'
#DATA=$1
#EVEN_FILES=${DATA}'/even/sep'${SEP}'/*I.uvGA'
#ODD_FILES=${DATA}'/odd/sep'${SEP}'/*I.uvGA'
#DIRNAME=$2
#EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
#ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
#CHAN='110_130'
#NBOOT=20
POL='I'
weight='L^-1'
WINDOW='none'
FRF='--frf'
NOFRFPATH='' #'--nofrfpath pspec128_uvGA/inject_sep'${SEP}'_0.01/pspec_pk_k3pk.npz' # path to one pspec_2d_to_1d.py output for NONFRF case
LMODE='' #'--lmode=12'
CHANGEC='--changeC'
#NGPS=5

### PSA64 Options ###

EVEN_FILES='/home/cacheng/capo/ctc/matt_data/even/*uvGAL'
ODD_FILES='/home/cacheng/capo/ctc/matt_data/odd/*uvGAL'
CALFILE='psa6240_v003'
CHAN='95_115'
SEP='0,1'
RA='.1_8.6'
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
DIRNAME='/home/saulkohn/pspec_banana_test_large/'
NBLG_MIN=1
NBLG_MAX=23
NLSTG_MIN=2
NLSTG_MAX=10
fi
#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}

# Stage 1: pspec_oqe_2d.py over range of injection levels
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,3,1)))"` ; do
    for NBLG in $(seq ${NBLG_MIN} ${NBLG_MAX}); do
        echo Making Directory ${DIRNAME}/nblg${NBLG}/
        mkdir ${DIRNAME}/nblg${NBLG}
        mkdir ${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}
        echo SIGNAL_LEVEL=${inject}
        
        # Stage 1: pspec_oqe_2d.py
        ${PATH2CAPO}/pspec_pipeline/pspec_oqe_2d.py ${LMODE} ${CHANGEC} --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} \
        -C ${CALFILE} -i ${inject} --weight=${weight} ${FRF} --output ${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject} \
        ${EVEN_FILES} ${ODD_FILES} --NGPS=${NBLG} 

        # Stage 2: pspec_2d_to_1d.py
        if [ ${NBLG} -eq ${NBLG_MIN} ]; then
            for NLSTG in $(seq ${NLSTG_MIN} ${NLSTG_MAX}); do
                echo Making Directory ${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/nlst${NLSTG} 
                mkdir ${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/nlst${NLSTG}
                ${PATH2CAPO}/pspec_pipeline/pspec_2d_to_1d.py --NGPS_LST=${NLSTG} \
                --output=${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/nlst${NLSTG}/ \
                ${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/pspec_oqe_2d.npz
            done
            ${PATH2CAPO}/pspec_pipeline/pspec_2d_to_1d.py --NGPS_LST=1 --output=${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/ \
            ${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/pspec_oqe_2d.npz
        else
            ${PATH2CAPO}/pspec_pipeline/pspec_2d_to_1d.py --NGPS_LST=1 --output=${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/ \
            ${DIRNAME}/nblg${NBLG}/inject_sep${SEP}_${inject}/pspec_oqe_2d.npz
        fi
    done
done

# Stage 3: Plot things
#${PATH2CAPO}/pspec_pipeline/nGroupPlots.py --path2data=${DIRNAME} --nblg_min=${NBLG_MIN} --nblg_max=${NBLG_MAX} --nlst_min=${NLSTG_MIN} --nlstmax=${NLSTG_MAX} --plotAll

