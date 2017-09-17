#!/bin/bash

[ "$USER" == "pnef" ]     && WorkDir=/u/at/pnef/Work/Code/Reclustering/
[ "$USER" == "swiatlow" ] && WorkDir=/u/at/swiatlow/nfs/projects/Reclustering/
[ "$USER" == "bpn7" ]     && WorkDir=/nfs/slac/g/atlas/u01/users/bnachman/SLAC_pythia/ForwardForldTest/
[ "$USER" == "jdamp" ]  && WorkDir= #Johannes fill in!
# add similar line if you are not pnef

SubFileLoc=`pwd`/_batchSingleSub.sh
#rm $SubFileLoc
DateSuffix=`date +%Y%m%d_%Hh%Mmin`

echo '#!/bin/bash
echo CD to $1
echo CMD is $2

cd $1
source setup.sh
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift
shift
echo Calling $cmd $*
$cmd $*
cp -r $JOBFILEDIR/*.root $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
#Process=2
pThatMin=160
pThatMax=320
BosonMass=800
#for Process in 1,2; do
for Process in `seq -4 -4`; do #1 3
    echo $Process
    for q in `seq 1 1`; do #15
	Queue=long
	nevents=1000
	njobs=100
	[ "$q" == "1" ]  && pThatMin=200 && pThatMax=400 
	[ "$q" == "2" ]  && pThatMin=1000 && pThatMax=1500
	[ "$q" == "3" ]  && pThatMin=1500 && pThatMax=2000
	[ "$q" == "4" ]  && pThatMin=2000 && pThatMax=2500
	[ "$q" == "5" ]  && pThatMin=2500 && pThatMax=3000
	[ "$q" == "6" ]  && pThatMin=3000 && pThatMax=3500
	LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_bsub_${mu}_
	OutDirFinal=`pwd`/files/${DateSuffix}
	mkdir -p `dirname $LogPrefix`
	mkdir -p $OutDirFinal
	echo
	echo "Submitting $njobs jobs each with $nevents events to $Queue"
	echo $LogPrefix
	for (( ii=1; ii<=$njobs; ii++ )) ;  do
            echo $ii
            OutDir=/scratch/${DateSuffix}_${ii}/
            bsub -q ${Queue} -R 'select[(!preempt&&rhel60&&cvmfs&&inet)]' -o $LogPrefix${ii}.log $SubFileLoc           \
		${WorkDir} ${OutDir} ${OutDirFinal} ./myexample.exe  \
		--OutFile ${OutDir}/Sample_mu_${mu}_nevents_${nevents}_job_${ii}_Process_${Process}_${pThatMin}_${pThatMax}.root \
		--NEvents ${nevents} \
		--Proc ${Process} \
		--pThatMin ${pThatMin} \
		--pThatMax ${pThatMax} \
	    
	done
    done
done
    
