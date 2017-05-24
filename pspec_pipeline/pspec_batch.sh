#!/bin/bash

### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>

#PSA64
if true
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
CALFILE='psa6622_v003'
RA='4_10'
SEP='0,2'
DATA=$1
EVEN_FILES=${DATA}'/even/sep'${SEP}'/*I.uvGA'
ODD_FILES=${DATA}'/odd/sep'${SEP}'/*I.uvGA'
DIRNAME=$2
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
CHAN='110_130'
#NBOOT=20
POL='I'
weight='L^-1'
WINDOW='none'
FRF='--frf'
NOFRFPATH='--nofrfpath pspec128_uvGA/inject_sep'${SEP}'_0.01/pspec_pk_k3pk.npz' # path to one pspec_2d_to_1d.py output for NONFRF case
LMODE='' #'--lmode=12'
CHANGEC='' #'--changeC'
fi
### PSA64 Options ###

#EVEN_FILES='/home/cacheng/capo/ctc/matt_data/lstbin_psa64_data_optimal/even/*uvGAL'
#ODD_FILES='/home/cacheng/capo/ctc/matt_data/lstbin_psa64_data_optimal/odd/*uvGAL'
#CALFILE='psa6240_v003'
#CHAN='95_115'
#SEP='0,1'
#RA='.1_8.6'
#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}

# Stage 1: pspec_oqe_2d.py over range of injection levels
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,3,1)))"` ; do
    mkdir ${DIRNAME}/inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}

    #calculate error bars across 5 groups (the SUBSET)
    ~/capo/pspec_pipeline/pspec_oqe_2d.py ${LMODE} ${CHANGEC} --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} \
    -C ${CALFILE} -i ${inject} --weight=${weight} ${FRF} --output ${DIRNAME}/inject_sep${SEP}_${inject} \
    ${EVEN_FILES} ${ODD_FILES} --NGPS=5

    #get the deepest power spectrum (the FULL)
    ~/capo/pspec_pipeline/pspec_oqe_2d.py ${LMODE} ${CHANGEC} --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} \
    -C ${CALFILE} -i ${inject} --weight=${weight} ${FRF} --output ${DIRNAME}/inject_sep${SEP}_${inject} \
    ${EVEN_FILES} ${ODD_FILES} --NGPS=1

    # Stage 2: pspec_2d_to_1d.py
    #(variance across the subset)
    ~/capo/pspec_pipeline/pspec_2d_to_1d.py ${FRF} ${NOFRFPATH} \
    --output=${DIRNAME}/inject_sep${SEP}_${inject}/ ${DIRNAME}/inject_sep${SEP}_${inject}/pspec_oqe_2d.npz
    
    #average over bls and lsts (todo: remove unnecessary bootstrap)
    ~/capo/pspec_pipeline/pspec_2d_to_1d.py ${FRF} ${NOFRFPATH} \
    --output=${DIRNAME}/inject_sep${SEP}_${inject}/ ${DIRNAME}/inject_sep${SEP}_${inject}/pspec_oqe_2d_final.npz

done
