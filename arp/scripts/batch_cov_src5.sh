#$ -S /bin/bash
ARGS=`pull_args.py $*`
SRCS=21:44:17.81_28:12:24.7,18:44:20.21_45:29:16.6,21:19:07.35_49:36:18.2,18:35:09.37_32:42:30.5,1:57:25.31_28:53:10.6,1:36:19.69_20:58:54.8,1:08:54.37_13:19:28.8,11:14:38.91_40:37:12.7,20:19:55.31_29:44:30.8,17:20:37.50_-0:58:11.6,0:42:53.96_52:07:34.8,18:56:36.10_1:20:34.8,4:08:02.40_43:00:31.1,10:01:31.41_28:48:04.0,9:21:18.65_45:41:07.2,15:04:55.31_26:01:38.9,5:42:50.23_49:53:49.1,21:55:53.91_37:55:17.9,16:28:35.62_39:32:51.3,21:23:54.38_25:02:07.2,14:11:21.08_52:07:34.8,1:37:22.97_33:09:10.4,8:13:17.32_48:14:20.5,22:45:49.22_39:38:39.8,20:14:17.81_23:33:42.8,4:18:02.81_38:00:58.6,16:51:05.63_5:00:17.4,3:19:41.25_41:30:38.7,5:04:48.28_38:06:39.8,4:37:01.87_29:44:30.8,vir,crab,cyg,cas,Sun
#SRCS=vir,crab,cyg,cas,Sun

echo cov_src5.py -C pgb015_v005 -s $SRCS -c 240_720_4 -x 4 -a cross,-1,-7,-13,-15 -p yy -r 5 -d 5 --clean=1e-2 --maxiter=100 $ARGS
cov_src5.py -C pgb015_v005 -s $SRCS -c 240_720_4 -x 4 -a cross,-1,-7,-13,-15 -p yy -r 5 -d 5 --clean=1e-2 --maxiter=100 $ARGS
