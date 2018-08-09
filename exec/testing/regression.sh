#!/bin/bash

usage () {
read -r -d '' help  <<-"EOT"
Perform some basic regression testing in a user's local repository.
IMPORTANT: Execute from decks directory only.

Usage: ../exec/testing/regression.sh [RUNSRC1 ... RUNSRCN -c --clean --help]

Where [...] are optional arguments representing rundeck template names
and additional flags:
  -c : perform compile-only verification
  --clean : clean all temporary files and directories created by this script
  --help : print this screen

Examples:

  Run without arguments:

  A)   ../exec/testing/regression.sh
  
  will regression test nonProduction_E_AR5_C12 and will write all results
  in decks directory.

  B)   ../exec/testing/regression.sh E4F40 E4TcadiF40 Earobio_g6c -c

  will compile-only E4F40, E4TcadiF40 and Earobio_g6c and all work will 
  be done in the decks directory.

  C)   ../exec/testing/regression.sh E4F40 E4TcadiF40 E4TctomasF40

  will run restart-regression on E4F40, E4TcadiF40 and E4TctomasF40.

  D) ../exec/testing/regression.sh --clean

  will remove ALL the temporary files and directories created by the
  script.

Errors, if any, will be printed on STDOUT.

Caveats:

A) User must preload working env modules and set MODELERC
   - Builds with compiler specified in MODELERC
   - Builds without optimization ("-O0 -g")
B) Runs interactively
   - OK for small rundecks and quick "sanity checks"
C) No separate scratch space
   - Builds and runs proceed in the decks directory
   - Beware of quotas
   - Builds and runs proceed sequentially
D) No baseline testing is performed - just internal consistency
   - serial vs 4 pes
   - checkpoint/restart vs continuous
E) Needs python version 2.7.x

EOT
   echo -e "$help"
   exit 1
}

clean () {
   local dirs=()
   rm -rf *.mk *.R *.diff *.log *cfg* *_bin templ
   dirs=`find . -maxdepth 1  -type l -exec ls -d {} \;`
   if [ ! -z "${dirs}" ]; then
      for d in "${dirs[@]}"; do
         rm -rf $(readlink -q ${d})
      done
   fi
   find . -type l -exec rm {} \;
   exit 1
}

root=`pwd`
scripts=$root/../exec/testing

if [ -z "$MODELERC" ]; then
   echo "Please set MODELERC."
   exit 1
else
   compiler=`grep COMPILER $MODELERC | awk -F= '{print $2}'`
fi

cnt=0
verification=restartRun
if [ "$#" -gt 0 ]; then
   if [ "$1" == "--clean" ]; then
      clean
   fi
   if [ "$1" == "--help" ]; then
      usage
   fi
   userArgs=( "$@" )
   rundecks=()
   for arg in "${userArgs[@]}"; do
      if [ "$arg" == "-c" ]; then
	 verification=compileOnly
      else
	 rundecks=( "${rundecks[@]}" "$arg" )
      fi
   done
else
   rundecks=("nonProduction_E_AR5_C12")
fi

   node=`uname -n`
   # We need the right python version on DISCOVER
   if [[ "$node" =~ discover || "$node" =~ dali || "$node" =~ borg ]]; then
      export PATH=/usr/local/other/SSSO_Ana-PyD/2.1.0/bin:$PATH
   fi

   repo=${root%/*}
   for run in "${rundecks[@]}"; do
      cp  $scripts/template.cfg $run.cfg
      sed -i -e "s|COMPILER|${compiler}|g" $run.cfg
      sed -i -e "s|REPO|${repo}|g" $run.cfg
      sed -i -e "s|MODELERC|${MODELERC}|g" $run.cfg
      sed -i -e "s/RUNDECK/${run}/g" $run.cfg
      sed -i -e "s/VERIFICATION/${verification}/g" $run.cfg
   done

python $scripts/regression.py ${rundecks[@]}
