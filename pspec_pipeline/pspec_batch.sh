#!/bin/bash

### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>


#### Carina's PSA128 Options ###
if True
then
DATA=$1
DIRNAME=$2
RA='1_10'
CALFILE='psa6622_v003'
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
SEP='0,2'
CHAN='110_130'
NBOOT=20
POL='I'
weight='I'
WINDOW='none'
FRF='--frf'
LMODE=''
CHANGEC=''
else
### psa64 Paths - enterprise ###
DATA=$1
EVEN_FILES=${DATA}'/even/sep0,1/*.uvGA'
ODD_FILES=${DATA}'/odd/sep0,1/*.uvGA'
DIRNAME=$2

## psa64 Options ###
RA='-0.1_8.6'
CALFILE='psa6240_v003'
SEP='0,1 1,1 -1,1'
CHAN='95_115'
NBOOT=60
POL='I'
weight='I'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='' #'--changeC'
fi

### Carina's PSA64 Options ###

#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}

for sep in $SEP; do
    mkdir ${DIRNAME}/sep${sep}
    EVEN_FILES=${DATA}'/even/sep'${sep}'/*.uvGAL'
    ODD_FILES=${DATA}'/odd/sep'${sep}'/*.uvGAL'
    EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
    ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
    # Stage 1: pspec_oqe_2d.py over range of injection levels
    for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,4,60)))"` ; do
        mkdir ${DIRNAME}/sep${sep}/inject_sep${sep}_${inject}
        echo SIGNAL_LEVEL=${inject}
        ~/src/capo/pspec_pipeline/pspec_oqe_2d.py --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} -C ${CALFILE} -b ${NBOOT} -i ${inject} ${dblfrf} --weight=${weight} --output ${DIRNAME}/sep${sep}/inject_sep${sep}_${inject} ${EVEN_FILES} ${ODD_FILES}
    
        # Stage 2: pspec_2d_to_1d.py
        ~/src/capo/pspec_pipeline/pspec_2d_to_1d.py --output ${DIRNAME}/sep${sep}/inject_sep${sep}_${inject}/ ${DIRNAME}/sep${sep}/inject_sep${sep}_${inject}/*boot*
    done
done

for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,4,60)))"` ; do
    mkdir ${DIRNAME}/inject_${inject}
    python ~/src/capo/pspec_pipeline/pspec_average_seps.py --outfile ${DIRNAME}/inject_${inject}/pspec_pk_k3pk.npz ${DIRNAME}/sep*/inject_sep*_${inject}/pspec_pk_k3pk.npz
done

python ~/src/capo/pspec_pipeline/pspec_final_confidence.py ${DIRNAME}/inject_*/pspec_pk_k3pk.npz --outfile ${DIRNAME}
