#!/bin/sh

cd $HOME/workspace/SnpSift/

mvn assembly:assembly
cp target/snpSift-1.6-jar-with-dependencies.jar $HOME/snpEff/SnpSift.jar

cd - 


