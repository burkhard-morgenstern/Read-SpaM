
all: 
		g++ -O3 -std=c++11 Bucket.cpp Seed.cpp readspam.cpp pattern.cpp Sequence.cpp Word.cpp -o readspam -fopenmp
 
