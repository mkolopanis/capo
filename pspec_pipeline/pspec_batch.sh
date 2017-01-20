#!/bin/bash
### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>#   ex: pspec_batch.sh /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/lstbin pspec_jan2017

### My Paths ###
DATA=$1
EVEN_FILES=${DATA}'/even/sep0,2/*I.uvGAL'
ODD_FILES=${DATA}'/odd/sep0,2/*I.uvGAL'
DIRNAME=$2

### My Options ###
RA='4_10'
CALFILE='psa6622_v003'
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
SEP='0,2'
CHAN='79_99'
NBOOT=20
POL='I'
weight='L^-1'
WINDOW='none'
FRFEOR='--frfeor' #to FRF the injected EOR, leave this on
SUBPCV='--sub_pCv' #to subtract pCv before bootstrapping, leave this on

#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}

# Stage 1: pspec_oqe_2d.py over range of injection levels
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,3,10)))"` ; do
    mkdir ${DIRNAME}/inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
    ~/src/capo/pspec_pipeline/pspec_oqe_2d.py --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} -C ${CALFILE} -b ${NBOOT} ${FRFEOR} -i ${inject} --weight=${weight} --output ${DIRNAME}/inject_sep${SEP}_${inject} ${EVEN_FILES} ${ODD_FILES}

    # Stage 2: pspec_2d_to_1d.py
    ~/src/capo/pspec_pipeline/pspec_2d_to_1d.py ${SUBPCV} --output ${DIRNAME}/inject_sep${SEP}_${inject}/ ${DIRNAME}/inject_sep${SEP}_${inject}/*boot*
done
