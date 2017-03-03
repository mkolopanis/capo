#$ -S /bin/bash
#$ -V
#$ -cwd 
#$ -l h_vmem=16G
#$ -l paper
#$ -o grid_output
#$ -e grid_output
#$ -N VIS_SIM

myargs=`pull_args.py $*`

echo my times: ${myargs}

name=`echo ${myargs} | cut -d " " -f 1`
echo first arg: ${name}

vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.1 --nchan 203 --inttime 31.65 --map gsm --mappath /home/cacheng/capo/ctc/images/gsm/gsm203/ --filename gsm_${name}.uv -C psa6622_v003 -a 1_4,1_47,1_68,1_19,1_75,1_48,1_18 ${myargs}

echo vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.1 --nchan 203 --inttime 31.65 --map gsm --mappath /home/cacheng/capo/ctc/images/gsm/gsm203/ --filename gsm_${name}.uv -C psa6622_v003 -a 1_4,1_47,1_68,1_19,1_75,1_48,1_18 ${myargs}



#startjd=2454500
#endjd=2454501
#inttime=20000
#numchunks=4 #number of processors... MUST match qsub command

#ARGS=`python -c "import numpy; import aipy; print ' '.join(map(str,numpy.arange(${startjd},${endjd},${inttime}/aipy.const.s_per_day)))"`
#echo times: ${ARGS}

#numtimes=`python -c "import numpy; import aipy; print len(numpy.arange(${startjd},${endjd},${inttime}/aipy.const.s_per_day))"`
#echo numtimes: ${numtimes}

#chunksize=$((${numtimes}/${numchunks}+1))
#echo numchunks: ${numchunks}
#echo chunksize: ${chunksize}

#timeindex=`python -c "import numpy; print ' '.join(map(str,numpy.arange(1,${numtimes}+1,${chunksize})))"`

#loopargs=`pull_args.py ${timeindex}`

#for chunk in ${loopargs}; do
    
    #echo ${chunk}
    #list="$(echo ${ARGS} | cut -d " " -f ${chunk}-$((${chunk}+${chunksize}-1)))" #list of times broken into chunks
    #echo mylist: ${list}
    #name=`echo ${list} | cut -d " " -f 1` 
    #echo ${name}
    #echo vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.1 --nchan 203 --inttime ${inttime} --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec203lmax200/ --filename /data2/home/cacheng/capo/ctc/tables/203files/pspec_${name}.uv -C psa898_v003 -a 0_16 ${list}
    #vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.1 --nchan 203 --inttime ${inttime} --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec203lmax200/ --filename /data2/home/cacheng/capo/ctc/tables/203files/pspec_${name}.uv -C psa898_v003 -a 0_16 ${list}
    #vis_simulation_v4.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime ${inttime} --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec100lmax100/ --filename test_${name}.uv -C psa898_v003 -a 0_16 ${list}
#done


