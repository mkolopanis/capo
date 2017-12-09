#!/bin/bash

### Sample Call ###
#   $ pspec_batch.sh <path to LST-binned files> <directory name to save all outputs>

echo "Welcome to the power spectrum pipeline!"

#PSA64
if true
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
CHAN='30_146 78_156'
#CHAN=' 30_50 51_71 78_98 95_115 103_123 127_147' #psa64 multiz bands
NBOOT=20
POL='I'
weight='I'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='--changeC' #throw out off diagonal terms of covariance.
NEIG='--neig=0'
NGPS=5
NGPS_LST=2
VERSION=2
#RMBLS='0_7,3_49,4_5,5_14,6_42,6_62,7_40,11_15,12_14,13_24,16_33,16_36,17_48,18_40,21_43,42_61,44_52,45_53,60_62'
#RMBLS='0_46,1_3,1_18,2_15,2_20,4_25,8_21,9_61,10_58,12_26,13_30,17_32,22_63,23_56,25_41,27_28,27_59,29_55,32_54,33_63,34_57,35_58,38_54,39_44,47_48,55_56'

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
weight='I' #'L^-1'
WINDOW='none'
FRF='--frf'
LMODE='' #'--lmode=12'
CHANGEC='' #'--changeC'
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

for chan in ${CHAN}; do
    if [ $(wc -w <<< ${CHAN}) -gt 1 ]; then
        mkdir -p ${DIRNAME}/$chan
        out_name=${chan}/inject_sep${SEP}
    else
        out_name=chan_${chan}_inject_sep${SEP}
    fi
    
    # Stage 1: pspec_oqe_2d.py over range of injection levels
    for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-3,5,1)))"` ; do
        out_dir=${out_name}_${inject}
        mkdir -p ${DIRNAME}/${out_dir}
        echo SIGNAL_LEVEL=${inject}
    
        ~/src/capo/pspec_pipeline/pspec_oqe_2d.py ${LMODE} ${CHANGEC} --window=${WINDOW} -a cross -p ${POL} -c ${chan} ${NEIG}\
        -C ${CALFILE} -i ${inject} --weight=${weight} ${FRF} --output ${DIRNAME}/${out_dir} -b ${NBOOT} \
        ${EVEN_FILES} ${ODD_FILES} --NGPS=${NGPS} --rmbls=${RMBLS}
    
        # Stage 2: pspec_2d_to_1d.py
        ~/src/capo/pspec_pipeline/pspec_2d_to_1d.py \
        --output=${DIRNAME}/${out_dir}/ --NGPS_LST=${NGPS_LST} -v ${VERSION} ${DIRNAME}/${out_dir}/*bootsigloss*.npz
        
    done
done
