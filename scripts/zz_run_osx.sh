#!/bin/sh
echo on
DIR="$( cd "$( dirname "$0" )" && pwd )"
cd ${DIR}
cp ../../_build/tcsam02 ./tcsam02
./tcsam02  -rs -nox  -configFile ../MCI.inp -phase 5  -calcOFL  -pin ../tcsam02.pin -nohess   

