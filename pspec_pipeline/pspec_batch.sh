#!/bin/bash

### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>

#### Carina's PSA128 Options ###
DATA=$1
EVEN_FILES=${DATA}'/even/sep0,2/*I.uvGAL'
ODD_FILES=${DATA}'/odd/sep0,2/*I.uvGAL'
DIRNAME=$2
RA='1_10'
CALFILE='psa6622_v003'
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
SEP='0,2'
CHAN='110_130'
NBOOT=20
POL='I'
weight='L^-1'
WINDOW='none'
FRF='--frf'
LMODE=''
CHANGEC=''

### psa64 Paths - enterprise ###
#DATA=$1
#EVEN_FILES=${DATA}'/even/sep0,1/*.uvGA'
#ODD_FILES=${DATA}'/odd/sep0,1/*.uvGA'
#DIRNAME=$2

### psa64 Options ###
#RA='-0.1_8.6'
#CALFILE='psa6240_v003'
#EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
#ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
#SEP='0,1'
#CHAN='95_115'
#NBOOT=20
#POL='I'
#weight='I'
#WINDOW='none'
#FRF='--frf'
#LMODE='' #'--lmode=12'
#CHANGEC='' #'--changeC'

### Carina's PSA64 Options ###

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
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,3,10)))"` ; do
    mkdir ${DIRNAME}/inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
    ~/capo/pspec_pipeline/pspec_oqe_2d.py ${LMODE} ${CHANGEC} --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} -C ${CALFILE} -b ${NBOOT} -i ${inject} --weight=${weight} ${FRF} --output ${DIRNAME}/inject_sep${SEP}_${inject} ${EVEN_FILES} ${ODD_FILES}

    # Stage 2: pspec_2d_to_1d.py
    ~/capo/pspec_pipeline/pspec_2d_to_1d.py --output ${DIRNAME}/inject_sep${SEP}_${inject}/ ${DIRNAME}/inject_sep${SEP}_${inject}/*boot*
done
#~/src/capo/pspec_pipeline/pspec_final_confidence.py ${DIRNAME}/inject_sep*/*npz
