#!/bin/bash

### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>

echo "Welcome to the power spectrum pipeline!"

#PSA64
if false
then
echo "Danny and Matt PSA64!"
CALFILE='psa6240_v003'
RA='.5_8.6'
SEP='0,1'
DATA=$1
EVEN_FILES=${DATA}'/even/sep'${SEP}'/*.uvGAL'
ODD_FILES=${DATA}'/odd/sep'${SEP}'/*.uvGAL'
DIRNAME=$2
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
CHAN='30_50'
#CHAN=' 30_50 51_71 78_98 95_115 103_123 127_147' #psa64 multiz bands
NBOOT=20
POL='I'
weight='L^-1'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='--changeC' #throw out off diagonal terms of covariance.
NGPS=5
NGPS_LST=2
VERSION=2

else
if false
then
echo "CARINA PSA128!"
### PSA128 options ###
CALFILE='psa6622_v003'
RA='4_10'
SEP='0,2'
DATA=$1
EVEN_FILES=${DATA}'/even/sep'${SEP}'/*I.uvGAL'
ODD_FILES=${DATA}'/odd/sep'${SEP}'/*I.uvGAL'
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
NBOOT=20
NGPS=5
NGPS_LST=2
VERSION=2

else
echo "CARINA PSA64!"
### PSA64 Options ###
POL='I'
weight='I' #'I' #'L^-1'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='--changeC'
NBOOT=20 # use 1 if doing version 4 (pspec_banana)
NGPS=5
NGPS_LST=2 # only matters for version 4 (otherwise it's not used)
VERSION=2 # version 4 is pspec_banana
EVEN_FILES='/home/cacheng/capo/ctc/matt_data/even/*uvGAL'
ODD_FILES='/home/cacheng/capo/ctc/matt_data/odd/*uvGAL'
CALFILE='psa6240_v003'
CHAN='95_115'
SEP='0,1'
RA='0.5_8.6'
RMBLS=''
EVEN_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C ${CALFILE} --ra=${RA} ${ODD_FILES[@]}`
DIRNAME=$2
fi
fi
#-----------------------------------------------------------------

# Make Power Spectrum Directory
mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}

for mode_num in `python -c "import numpy; print ' '.join(map(str,numpy.arange(0,22,1)))"` ; do
    mkdir -p ${DIRNAME}/project_${mode_num}_modes
    out_dir=project_${mode_num}_modes/inject_sep${SEP}
    echo MODE_NUM=${mode_num}
 
    # Stage 1: pspec_oqe_2d.py over range of injection levels
    for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-3,3,20))) + ' ' + ' '.join(map(str, -numpy.logspace(-3,3,20)))"` ; do
        out_dir_pspec=${out_dir}_${inject}
        mkdir -p ${DIRNAME}/${out_dir_pspec}
        echo SIGNAL_LEVEL=${inject}
    
        ~/capo/pspec_pipeline/pspec_oqe_2d.py ${LMODE} ${CHANGEC} --window=${WINDOW} -a cross -p ${POL} -c ${CHAN} -C ${CALFILE} -i ${inject} --weight=${weight} ${FRF} --output ${DIRNAME}/${out_dir_pspec} -b ${NBOOT} --NGPS=${NGPS} --rmbls=${RMBLS} --mode_num=${mode_num} ${EVEN_FILES} ${ODD_FILES}
    
        # Stage 2: pspec_2d_to_1d.py
        ~/capo/pspec_pipeline/pspec_2d_to_1d.py \
        --output=${DIRNAME}/${out_dir_pspec}/ --NGPS_LST=${NGPS_LST} -v ${VERSION} ${DIRNAME}/${out_dir_pspec}/*bootsigloss*.npz
       
    done
        
    # Stage 3: pspec_final_sigloss_dist_v2.py
    cd ${DIRNAME}/project_${mode_num}_modes
    ~/capo/pspec_pipeline/pspec_final_sigloss_dist_v2.py --sep=${SEP}
    cd ../..

done
