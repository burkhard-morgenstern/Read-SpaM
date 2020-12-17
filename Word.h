#ifndef WORD_H_
#define WORD_H_

#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream> 

class Word{
    public:
        Word() = default;
		Word(uint32_t, uint32_t);

        uint32_t key;
        uint32_t pos;

		uint32_t getKey();
		uint32_t getPos();
		uint32_t* getPosPointer();

		void setKey(uint32_t);

		bool operator<( const Word& val ) const { 
			return key < val.key; 
		}
};
#endif