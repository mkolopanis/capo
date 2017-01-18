### My Data ###
PREFIX='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/lstbin/'
EVEN_FILES=${PREFIX}'even/sep0,2/*I.uvGAL'
ODD_FILES=${PREFIX}'odd/sep0,2/*I.uvGAL'

### My Options ###
RA='4_10' #'1_10'
EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C psa6622_v003 --ra=${RA} ${ODD_FILES[@]}`
SEP='0,2'
CALFILE='psa6622_v003'
CHAN='79_99'
NBOOT=20
POL='I'
WINDOW='none'
DIRNAME='TEST_pspec_'${SEP}'_'${CHAN}'_'${RA}'_'${POL}
FRFEOR='--frfeor' #to FRF the injected EOR, leave this on

#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}
cd ${DIRNAME}

# Stage 1: pspec_oqe_2d.py over range of injection levels
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-1,3,10)))"` ; do
    mkdir inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
    ~/capo/pspec_pipeline/pspec_oqe_2d.py --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} -C ${CALFILE} -b ${NBOOT} ${FRFEOR} -i ${inject} ${EVEN_FILES} ${ODD_FILES}
    mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
done

# Stage 2: pspec_2d_to_1d.py 
~/capo/pspec_pipeline/pspec_2d_to_1d.py 
