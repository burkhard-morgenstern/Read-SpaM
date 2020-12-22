#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream> 
#include <algorithm> 
#include <omp.h>
#include <map>
#include <libgen.h>
#include "Word.h"
#include "Seed.h"
#include "Bucket.h"

class Sequence{

	private:
		uint32_t wordCount;
		std::vector<char> seq;
		std::vector<char> seqRev;
		std::vector<Word> sortedWords;
		std::vector<Word> sortedWordsRev;
		std::vector<uint32_t> firstBuckets;
		std::vector<uint32_t> firstBucketsRev;
		std::string header;

		void cleanup();
	public:

		Sequence(std::string, std::string &);

		std::string getHeader(){
			return header;
		}

		std::vector<char>& getSequence(){
			return seq;
		}
		void sortFirstBits(Seed &);
		void sortNextBits(Seed &);

		void sortFirstBitsRev(Seed &);
		void sortNextBitsRev(Seed &);

		/*
		 * Reads in first bits of words while sampling using the
		 * specified ratio. 'reverse' specifies the strand.
		 */
        void sortFirstBits(Seed& seed, double ratio, bool reverse, int seedMinhash);

		double compareSequences(Sequence & qry, Seed & s, int threads, int threshold, bool writeHistogram, std::string outputFolder);

		static std::vector<Sequence> read(std::string);
};

struct MatchScore{ 
	uint32_t posS1;
	uint32_t posS2;
	int32_t score;
	char mismatches;

	bool operator<( const MatchScore& val ) const { 
		return score > val.score; 
	}

	MatchScore(uint32_t p1, uint32_t p2, int32_t s, char mm)
	{
		posS1 = p1;
		posS2 = p2;
		score = s;
		mismatches = mm;
	}
};

#endif