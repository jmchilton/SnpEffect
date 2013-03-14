

									Release process
									---------------

1) Run JUnit tests (TestSuiteAll) and make sure all of them pass

2) Run integration tests

3) Change version numbers:

	- Update SnpEff pom.xml 
	- Update SnpSift pom.xml
	- Update scripts_build/config.sh 

4) Build JAR files, download databases, build databases, etc.
	
	cd ~/snpEff
	./scripts_build/build.sh 

5) Update Galaxy's snpEff.xml
		# Create galaxy menu (add output to galaxy/snpEff.xml)
		java -jar snpEff.jar cfg2table galaxy 

6) Upload to Galaxy ToolShed: http://toolshed.g2.bx.psu.edu/
		Reference: http://wiki.g2.bx.psu.edu/Tool%20Shed

-------------------------------------------------------------------------------		
		
How to add the libraries to your local Maven repository:

    http://maven.apache.org/general.html#importing-jars

    mvn install:install-file -Dfile=picard-1.77.jar -DgroupId=net.sf.picard -DartifactId=Picard -Dversion=1.77 -Dpackaging=jar
    mvn install:install-file -Dfile=sam-1.77.jar -DgroupId=net.sf.samtools -DartifactId=Sam -Dversion=1.77 -Dpackaging=jar
   
