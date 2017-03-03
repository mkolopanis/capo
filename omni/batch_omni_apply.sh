#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l scratch=9G
#$ -l paper
#$ -N OMNI_APPLY
ARGS=`pull_args.py $*`

echo ${ARGS}

for f in ${ARGS}; do
    ~/capo/omni/omni_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v3_xtalk/%s.npz -p yy --xtalk --ba 7,16,27,50,54,56,103 --ubls="0,1;1,1;1,-1;0,2;0,3;1,2;1,-2" ${f} #S2E4yy
    #~/capo/omni/omni_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v3_xtalk/%s.npz -p xx --xtalk --ba 8,16,17,24,28,29,34,38,68,85 --ubls="0,1;1,1;1,-1;0,2;0,3;1,2;1,-2" ${f} #S2E4xx
    #~/capo/omni/omni_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v3_xtalk/%s.npz -p xx --xtalk --ba 8,24,34,38,85,107 --ubls="0,1;1,1;1,-1;0,2;0,3;1,2;1,-2" ${f} #S2E3xx
    #~/capo/omni/omni_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v3_xtalk/%s.npz -p yy --xtalk --ba 7,16,34,56 --ubls="0,1;1,1;1,-1;0,2;0,3;1,2;1,-2" ${f} #S2E3yy
    #~/capo/omni/omni_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/%s.npz -p xx --xtalk --ba 8,16,24,34,38,53,63,74,85 --ubls="0,1;1,1;1,-1;0,2;0,3;1,2;1,-2" ${f} #S2E2xx
    #~/capo/omni/omni_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/%s.npz -p yy --xtalk --ba 7,16,34,56 --ubls="0,1;1,1;1,-1;0,2;0,3;1,2;1,-2" ${f} #S2E2yy
    #if ((${f:33:7} < 2456679 )); then
    #    echo working on ${f}, which is in Epoch 1...
    #    ~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v3_xtalk/%s.npz -p yy --xtalk --ubls="0,2;0,1" --ba 2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,100,7,56,84 ${f} #S1E1yy
    #fi
    
    #if (( ${f:33:7} > 2456678 )); then
    #    echo working on ${f}, which is in Epoch 2...
    #    ~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/%s.npz -p yy --xtalk --ubls="0,2;0,1" --ba 100,7,56,84 ${f} #S1E2yy
    #    ~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk_removedegenOFF/%s.npz -p xx --xtalk --ubls="0,2;0,1" --ba 34,84,100 ${f} #S1E2xx
    #fi
done
