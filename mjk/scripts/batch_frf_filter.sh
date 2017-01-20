#!/bin/bash

seps='sep0,1 sep1,1 sep-1,1'
cal=psa6240_v003
appelation='uvGA'
chan='101'

indir=$1
outdir=$2

scriptsdir=/home/mkolopanis/src/capo/mjk/scripts

alietal=false

declare -A ants
ants[sep0,1]=0_44
ants[sep1,1]=1_3
ants[sep-1,1]=1_48

days='even odd'

printf 'This is Batch Fringe Rate Filter using:\n'
printf 'calfile: %s \n' $cal
sleep 1.5

if [ ! -d $outdir   ]; then
    printf 'This output director appears to be new\n'
    printf 'Creating Directoires now\n'
    for path in $days; do
        for sep in $seps; do
            mkdir -p ${outdir}/${path}/${sep}
        done
   done
fi

for path in $days; do
    for sep in $seps; do
        printf 'Checking for Old FRF Data\n'
        outpath=$outdir/$path/$sep
        if [[ $(ls -d $outpath/*.${appelation}L) ]] ; then
            printf 'Found %s folders\n' $(ls -d $outpath/*.${appelation}L| wc -l)
            printf 'First deleting old FRF Filtered Data \n'
            sleep 1
            printf 'I wil delete: \n'
            sleep 1
            printf '%s \n' $(ls -d $outpath/*.${appelation}L)
            read -p "Are you sure you want to delete these? [y/N] " -n 1 -r
            printf '\n'
            if [[ $REPLY =~ ^[Yy]$ ]]
                then
                    rm -rf $(ls -d $outpath/*.${appelation}L)
                    printf 'Files are deleted\nContinuting to Fringe Rate Filter\n'
                else
                    printf 'Nothing is deleted\ncontinuing\n'
                    continue
            fi
        fi

        files=$(ls -d $indir/$path/$sep/*.${appelation})
        printf 'Filtering %s by selecting ant %s \n' $sep ${ants[$sep]}

        if ${alietal}; then

            printf '%s/frf_filter.py -C %s  --alietal -a %s -c %s  --outpath=%s/' \
            $scriptsdir  $cal ${ants[$sep]} $chan h ${outpath}
            "${scriptsdir}/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal --alietal \
            $files --outpath=${outpath} -c $chan

        else
            printf '%s/frf_filter.py -C %s   -a %s -c %s --outpath=%s/' $scriptsdir $cal ${ants[$sep]} $chan ${outpath}
            "${scriptsdir}/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal\
            $files --outpath=${outpath} -c $chan --pol=I

        fi
    done
done
