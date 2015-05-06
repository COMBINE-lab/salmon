/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#ifndef __BIASINDEX_HPP__
#define __BIASINDEX_HPP__

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>


#include <algorithm>

template <class T>
void endswap(T *objp)
{
  unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
  std::reverse(memp, memp + sizeof(T));
}

template <class T>
void endswap(std::vector<T>& v) {
	for ( size_t i = 0; i < v.size(); ++i ) { endswap(&v[i]); }
}

class Indexer {
    public:
	Indexer(std::vector<double>::iterator it) : it_(it), hasData_(true) {};
	Indexer(double val) : val_(val), hasData_(false) {};

	inline double operator[](size_t offset) {
		return hasData_ ? *(it_+offset) : val_;
	}

    private:
    	std::vector<double>::iterator it_;
    	double val_;
    	bool hasData_;
};

class BiasIndex {

	using Bias = double;
	using Offset = size_t;

public:
	BiasIndex() : haveBias_(false) {}
	
	BiasIndex(const std::string base) : haveBias_(true) {
		auto dictName = base+".dict";
		auto binName = base+".bin";

		std::cerr << "reading bias dictionary ... ";
		using std::string;
		std::ifstream din(dictName, std::ios::in);
		
		size_t numBiasTerms;
		din >> numBiasTerms;

		string tname;
		Offset length;
		Offset offset=0;
		while ( din >> tname >> length ) {
			offsetMap_[tname] = offset;
			offset += length;
		}
		din.close();
		std::cerr << "done\n";

		std::cerr << "reading bias terms ... ";
		biases_.resize(numBiasTerms);

		std::ifstream bin(binName, std::ios::binary);
		bin.read( reinterpret_cast<char*>(&biases_.front()), sizeof(double)*numBiasTerms );
		bin.close();
		std::cerr << "done\n";
	}

	inline Indexer getBiases( const std::string& tname ) {
		return haveBias_ ? Indexer( biases_.begin() + offsetMap_[tname] ) : Indexer(1.0);
	}

private:
	bool haveBias_;
	std::vector<Bias> biases_;
	std::unordered_map<std::string, size_t> offsetMap_;
};

#endif //__BIASINDEX_HPP__