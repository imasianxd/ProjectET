# ProjectET
epitope prediction from a FASTA sequence using a combination of netMHCII-2.2, netMHCIIpan-3.0 and iedb software 

Advanced Informations:

Dependencies:
	mono-complete
	(can be installed using apt-get install mono-complete)
	EPPlus.dll
	Jia.dll
	Microsoft.Office.Interop.Excel.dll
	Newtonsoft.Json.dll
	OxyPlot.dll
	Renci.SshNet.dll
	(All DLL libraries can be found in ./projectET/dlls)

***you should never need to follow these steps unless instructed to do so***
To change the source code, open ./projectET/ProjectET/findEpitope.cs and make changes.
Then, you must build the program to apply the changes, I suggest building in debugging mode (default).
To build the program, open a terminal and navigate to /home/j/Desktop/ProjectET_Source or wherever the .sln file is.
type xbuild into terminal
ignore all warning.
If error exists during compiling, make sure you have mono-complete installed first by typing sudo apt-get install mono-complete and you did not make any syntax errors while changing the code.

