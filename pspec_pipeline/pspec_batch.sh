#!/bin/bash

### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>

echo "Welcome to the power spectrum pipeline!"

#PSA64
if false
then
echo "Danny and Matt PSA64!"
CALFILE='psa6240_v003'
RA='.1_8.6'
SEP='0,1'
DATA=$1
EVEN_FILES=${DATA}'/even/sep'${SEP}'/*.uvGAL'
ODD_FILES=${DATA}'/odd/sep'${SEP}'/*.uvGAL'
DIRNAME=$2
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
CHAN='95_115'
#NBOOT=20
POL='I'
weight='I'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='' # --changeC' #throw out off diagonal terms of covariance.
NGPS=1
NGPS_LST=2

else
if false
then
echo "CARINA PSA128!"
### PSA128 options ###
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
POL='I'
weight='L^-1'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='--changeC'
NGPS=1
NGPS_LST=2

else
echo "CARINA PSA64!"
### PSA64 Options ###
POL='I'
weight='L^-1'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='--changeC'
NGPS=1
NGPS_LST=2
EVEN_FILES='/home/cacheng/capo/ctc/matt_data/even/*uvGAL'
ODD_FILES='/home/cacheng/capo/ctc/matt_data/odd/*uvGAL'
CALFILE='psa6240_v003'
CHAN='95_115'
SEP='0,1'
RA='.1_8.6'
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
DIRNAME=$2
fi
fi
#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}

# Stage 1: pspec_oqe_2d.py over range of injection levels
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,3,1)))"` ; do
    mkdir ${DIRNAME}/inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}

    ~/capo/pspec_pipeline/pspec_oqe_2d.py ${LMODE} ${CHANGEC} --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} \
    -C ${CALFILE} -i ${inject} --weight=${weight} ${FRF} --output ${DIRNAME}/inject_sep${SEP}_${inject} \
    ${EVEN_FILES} ${ODD_FILES} --NGPS=${NGPS}

    # Stage 2: pspec_2d_to_1d.py
    ~/capo/pspec_pipeline/pspec_2d_to_1d.py \
    --output=${DIRNAME}/inject_sep${SEP}_${inject}/ --NGPS_LST=${NGPS_LST} ${DIRNAME}/inject_sep${SEP}_${inject}/pspec_oqe_2d.npz

done
