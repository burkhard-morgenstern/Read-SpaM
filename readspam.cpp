/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <iostream> 
#include <stdlib.h>
#include <getopt.h>
#include <iomanip> 
#include "Sequence.h"
#include "Seed.h"
 
void printHelp(){
 std::string help = 
    "\nUsage: ./fswm [options] <filelist>"
    "\n"
   	"\nformat:"  
	"\n\t	<filelist>:"
	"\n\t A plain text file, specifying the relative path to each input dataset."
	"\n\t To create the filelist simply run:"
	"\n\t 	ls ./path/to/input/* > filelist	"
	"\n\t 	(assuming all your proteome files (*.fasta, *.faa, etc.) are stored in the input-folder.)"
	"\n\t Sequence must be in FASTA format. Each genome must be in its own FASTA file."
	"\n\t  There can be multiple header in each FASTA file"
    "\n\t >Gene1"
    "\n\t ATAGTAGATGAT.."
    "\n\t >Contig2"
    "\n\t ATAGTAGATGAT.."
    "\n\t >Contig3"
    "\n\t ATGATGATGATGATG.."
    "\n\t .."
    "\n\t "
    "\nOptions:"        
    "\n\t -h: print this help and exit"
    "\n\t -k <integer>: pattern weight (default 12)"
    "\n\t -l <integer>: don't care positions (default 100)"
    "\n\t -t <integer>: numer of threads (default: 10)"
    "\n\t -s <integer>: the minimum score of a spaced-word match to be considered homologous (default: 0)"
    "\n";
	std::cout << help << std::endl;
}

int weight = 12;
int dontCare = 100;
int threads = 10;
int threshold = 0;

void parseParameters(int argc, char *argv[]){
	int option_char;
	 while ((option_char = getopt (argc, argv, "l:k:t:hs:")) != -1){ 
		switch (option_char){ 
			case 's': 
				threshold = atoi (optarg); 
				break;
			case 'l': 
				dontCare = atoi (optarg); 
				break;
			case 'k': 
				weight = atoi (optarg); 
				if(weight<8 || weight > 16){
					std::cerr << "Weight (-k) must be between 8 and 16"<< std::endl;
					exit (EXIT_FAILURE);
				}
				break;
			case 't': 
				threads = atoi (optarg); 
				if(threads<1){
					std::cerr << "threads (-t) must be an integer larger than 0"<< std::endl;
					exit (EXIT_FAILURE);
				}
				break;
			case 'h': 
				printHelp();
				exit (EXIT_SUCCESS);
				break;
			case '?': 
				printHelp();		
				exit (EXIT_FAILURE);
      	}
	}
}

void writeDmat(std::vector<std::vector<double> > dmat, std::vector<Sequence>& sequences){
	std::cout << sequences.size() << std::endl;
	for (int i = 0; i < sequences.size(); i++) 
	{
		std::string name = sequences[i].getHeader();
		for(int k = 0; k < 10; k++){
			if(k >= name.length())
				std::cout << " ";
			else
				std::cout << name[k];
		}
     	for (int j = 0; j < sequences.size(); j++) 
     	{
			if (i > j) 
	    			std::cout << std::fixed <<std::setprecision(12) << dmat[i][j] << "  ";
			else if(j>i)
				std::cout << std::fixed<< std::setprecision(12) << dmat[j][i] << "  ";
			else
					std::cout << std::setprecision(12) << "0" << "  ";
     	}
      		std::cout << std::endl;
	}
}

int main(int argc, char *argv[]){
	if(argc < 2)
	{
		printHelp();		
		exit (EXIT_FAILURE);
	}
	parseParameters(argc,  argv);
	std::string list(argv[argc-1]);
	Seed seed(weight,dontCare);
	std::vector<Sequence> sequences = Sequence::read(list);
	std::vector<std::vector<double> >DMat(sequences.size(), std::vector<double>(sequences.size(),0));
	Seed::init();
	omp_set_dynamic(0);     
	omp_set_num_threads(threads);
	if(sequences.size() < 2){
		std::cerr << "there must be at least 2 sequences"<< std::endl;
		exit (EXIT_FAILURE);
	}
	#pragma omp parallel for schedule(runtime)
	for(int i = 0; i < sequences.size();i++)
	{
		sequences[i].sortFirstBits(seed);
		sequences[i].sortFirstBitsRev(seed);
	}

	for(int i = 0; i < sequences.size(); i++)
	{
		sequences[i].sortNextBits(seed);
		sequences[i].sortNextBitsRev(seed);
	}
	int number = (sequences.size()*(sequences.size()-1))/2;
	int cnt = 1;
	for(int i = 0; i < sequences.size(); i++)
	{
		for(int j = i + 1; j < sequences.size(); j++)
		{
			DMat[i][j] = sequences[i].compareSequences(sequences[j], seed, threads, threshold);
			DMat[j][i] = DMat[i][j];
		}
	}
	writeDmat(DMat, sequences);
}
