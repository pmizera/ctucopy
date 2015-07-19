/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2012 Petr Fousek, FEE CTU Prague

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

#ifndef TYPES_H_
#define TYPES_H_

#include <stdint.h>

typedef uint16_t uSint;
typedef uint32_t uLint;
typedef int16_t Sint;
typedef int32_t Lint;

using namespace std;

// structures definition file
template <class T> class Vec {
public:
	T* vec;
	int size;
	Vec(int s) : vec(new T[s]), size(s) {
		for (int i=0; i<size; i++)
			vec[i] = 0;
	}
	~Vec() { delete [] vec; }
	int get_size() const { return size; }
	T& operator[] (int i) {
		// debugging code, to be removed:
//		if (i>=size) {cerr << "TYPES: Index out of range (i="<<i<<",size="<<size<<"!\n";throw;}
		// end of dbg code
		return vec[i];
	}
	T operator[] (int i) const { return vec[i]; }
	void print(int start) {
		for (int i=start; i<start+size; i++)
			cout << vec[i%size] << "\t";
		cout << endl;
	}
	//	friend ostream &operator<< <T>(ostream&, Vec<T>&);
};

template<class T> ostream& operator<<(ostream& c, Vec<T>&item) {
	//	c << "size =  " << item.size << endl;
	for (int i=0; i<item.size; i++) {c << item.vec[i] << "\t";}	c << endl;
	return c;
};


template <class T> class Mat {
  int r,c;          // rows and columns
  T **data;

public:
	Mat(int r_, int c_) :  r(r_), c(c_) {
		data = new T*[r];
		for (int i=0; i<r; i++) data[i] = new T[c];
		for (int i=0; i<r; i++)
			for (int j=1; j<c; j++)
				data[i][j] = 0;
	}
	~Mat() {
		for (int i=0; i<r; i++) delete data[i];
		delete [] data;
	}
	int rows() const { return r; }
  	int cols() const { return c; }
    T& operator()(int i, int j) { return data[i][j]; }
    void print(int start) {
        for (int i=start; i<start+r; i++) {
            for (int j=0; j<c; j++)
                cout << data[i%r][j] << "\t";
            cout << endl;
        }
        cout << endl;
    }
};

template<class T> ostream& operator<<(ostream& c, Mat<T>&item) {
	//	c << "size =  " << item.size << endl;
	for (int i=0; i<item.r; i++)
		for (int j=0; j<item.c; j++)
			c << item.data[i][j] << "\t";
		c << endl;
	return c;
};



#endif /*TYPES_H_*/
