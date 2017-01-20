#!/bin/bash
### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>#   ex: pspec_batch.sh /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/lstbin pspec_jan2017

### My Paths ###
DATA=$1
EVEN_FILES=${DATA}'/even/sep0,1/lst.*.uvGAL'
ODD_FILES=${DATA}'/odd/sep0,1/lst.*.uvGAL'
DIRNAME=$2

### My Options ###
RA='0_8.5'
CALFILE='psa6240_v003'
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
SEP='0,1'
CHAN='95_115'
NBOOT=30
POL='I'
WINDOW='none'
FRFEOR='--frfeor' #to FRF the injected EOR, leave this on
SUBPCV='' #to subtract pCv before bootstrapping, leave this on

#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}

# Stage 1: pspec_oqe_2d.py over range of injection levels
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,3,10)))"` ; do
    mkdir ${DIRNAME}/inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
    ~/src/capo/pspec_pipeline/pspec_oqe_2d.py --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} -C ${CALFILE} -b ${NBOOT} ${FRFEOR} -i ${inject} --output ${DIRNAME}/inject_sep${SEP}_${inject} ${EVEN_FILES} ${ODD_FILES}

    # Stage 2: pspec_2d_to_1d.py
    ~/src/capo/pspec_pipeline/pspec_2d_to_1d.py ${SUBPCV} --output ${DIRNAME} ${DIRNAME}/inject_sep${SEP}_${inject}/*boot*
done
