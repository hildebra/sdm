/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand
email: Falk.Hildebrand@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "IO.h"
#include "Benchmark.h"
#include "ReadMerger.h"


int main(int argc, char* argv[])
{
    bool test = false;
    if (test) {
//        std::cout << "test" << std::endl;
        //Merge::testMerge();
//
//        std::string file1 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/r1_30.fastq";
//        std::string file2 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/r2_30.fastq";

//        std::string file1 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/r1_sub.fastq";
//        std::string file2 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/r2_sub.fastq";

//        std::string file1 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/fatih_sm.1.fq";
//        std::string file2 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/fatih_sm.2.fq";

//        std::string file1 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/1.fq";
//        std::string file2 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/2.fq";

        std::string file1 = "C:/Users/hildebra/OneDrive/science/data/test/dada2Seed/testMrgR1.txt";
        std::string file2 = "C:/Users/hildebra/OneDrive/science/data/test/dada2Seed/testMrgR2.txt";

//        std::string file1 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/test_cases.fq";
//        std::string file2 = "/usr/users/QIB_fr017/fritsche/Projects/sdm2/data/mergetest/test_cases2.fq";
//
        std::ifstream is1(file1.c_str());
        std::ifstream is2(file2.c_str());
//
//        const std::string revtest = "hallodri";
//        std::cout << revtest << std::endl;
//        Merge::reverseStringInPlace(const_cast<char *>(revtest.c_str()), revtest.length());
//        std::cout << revtest << std::endl;
//
//        std::cout << "test merge with reads" << std::endl;
		ReadMerger* RM = new ReadMerger();
        RM->testMergeWithReads(is1, is2);

        exit(0);
    }



    Benchmark sdm_benchmark("Time taken: ");
    sdm_benchmark.start();
    
	if (argc<3){
		//help_options,help_map,help_commands
		if (argc==2){
			if (string(argv[1])=="-help_commands"){
				printCmdsHelp();
			} else if (string(argv[1])=="-help_options"){
				printOptionHelp();
			} else if (string(argv[1])=="-help_map"){
				printMapHelp();
			}
			else if (string(argv[1]) == "-version" || string(argv[1]) == "-v"){
				printVersion();
			}
			exit(0);
		}
		general_help();
		exit(0);
	}
	cdbg("DEBUG mode\n");

	ini_DNAconstants();

	Announce_sdm();

	OptContainer cmdArgs;
	readCmdArgs(argc, argv, cmdArgs);
	cdbg("CmdArgs read\n");

	//rewrites header names
	rewriteNumbers(cmdArgs);
	cdbg("CmdArgs modified\n");

	//reads the sdm_options file
	cdbg("Setting up Filter\n");
	//shared_ptr<Filters> fil = make_shared<Filters>(&cmdArgs);
	
	Filters* fil = new Filters(cmdArgs);
	cdbg("filter setup\n");

	
	//bool bReads = fil->readMap();
	bool bReads = fil->readMap();
	if (!bReads) { cerr << "Failed to read Map.\n"; exit(3); }
	cdbg("map is read\n");

	fil->setcmdArgsFiles();
	cdbg("CmdArgs set in filter\n");

	clock_t tStart = clock();
	
	
	//main function
	separateByFile(fil, cmdArgs);
	//end main function	
 
	delete fil;
 

//	fprintf(stderr,"Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	sdm_benchmark.stop();
	sdm_benchmark.printResults(std::cerr);
	
	return 0;
}



