#!/bin/sh

cd $HOME/workspace/SnpEff/

mvn assembly:assembly
cp target/snpEff-2.1a-jar-with-dependencies.jar $HOME/snpEff/snpEff.jar

cd - 


