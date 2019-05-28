echo "Please Wait...Fetching tarball from server"
wget "bfjia.net/ftp/ProjectET.tar.gz";
tar -xvf ./ProjectET.tar.gz;
echo "Cleaning up"
rm -rf ./ProjectET.tar.gz;
rm -rf ./Results
rm -rf ./Old_Runs
rm -rf ./Protein_Sequence
mkdir Protein_Sequence
rm -rf ./Modified_Protein_Sequence
mkdir Modified_Protein_Sequence
rm -rf ./iedb/results
rm -rf ./netMHCII-2.2/results
rm -rf ./netMHCII-2.2/modified_results
rm -rf ./netMHCIIpan-3.0/results
rm -rf ./Debug_Settings.json
echo "{\"stepsCompleted\":-1,\"cutoff\":0.1,\"analysisTime\":0.0,\"isError\":false,\"nonPubSequence\":false}" > Debug_Settings.json
echo "ReplaceMeWithAccession>ReplaceMeWithName" > Protein_Accession_Dictionary.txt
rm -rf ./accession.txt
rm -rf netMHCIITemp
mkdir netMHCIITemp
echo Done.
