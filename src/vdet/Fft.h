/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2012 FEE CTU Prague

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

#ifndef DSP_FFT_H
#define DSP_FFT_H

#include <cassert>
#include "SmartPtr.h"
#include "Complex.h"

namespace DSP
{

/*
**  Sample Iterator
**  ~~~~~~~~~~~~~~~
**
**  Sample Iterator must have defined following operators
**
**  a) operator==, operator!=  // (in)equality test
**
**  b) operator*  // dereference operator (should return Complex or double)
**
**  c) operator++ // prefix and postfix
**/

class BaseFft
{
  public: // Typedefs
    typedef Complex*        iterator;
    typedef Complex const*  const_iterator;

  public:
    int Points( void ) const
    { return _nPoints; }

    double Intensity( int i ) const
    {
      assert( i >= 0 );
      assert( i < _nPoints );
      return _x[i].Mod() / _sqrtPoints;
    }

    double Mod( int i ) const
    {
      assert( i >= 0 );
      assert( i < _nPoints );
      return _x[i].Mod();
    }

    double Re( int i ) const
    { 
      assert( i >= 0 );
      assert( i < _nPoints );
      return _x[i].Re();
    }

    double Im( int i ) const
    { 
      assert( i >= 0 );
      assert( i < _nPoints );
      return _x[i].Im();
    }

    Complex const& operator[]( int i ) const
    { 
      assert( i >= 0 );
      assert( i < _nPoints );
      return _x[i];
    }

    const_iterator Begin( void ) const
    { return &_x[0]; 
    }

    const_iterator End( void ) const
    { return &_x[0] + _nPoints; 
    }

    // Member template
    template<class SI>
    void Process( SI begin, SI end )
    {
      // Initialize the FFT buffer
      SI it = begin;
      int i;
      for ( i = 0; it != end && i < _nPoints; ++it, ++i )
      { _x[_bitRev[i]] = Complex(*it);
      }

      assert( it == end && i == _nPoints );

      // step = 2 ^ (level-1)
      // increm = 2 ^ level;
      int step = 1;
      for ( int level = 1; level <= _logPoints; ++level )
      {
        int increm = step << 1;
        for ( int j = 0; j < step; ++j )
        {
          // U = exp ( sign 2 PI j / 2 ^ level )
          Complex U = _w[level][j];
          for ( int i = j; i < _nPoints; i += increm )
          {
            // butterfly
            Complex T = U;
            T *= _x[i + step];
            _x[i + step] = _x[i];
            _x[i + step] -= T;
            _x[i] += T;
          }
        }

        step <<= 1;
      }
    }

  protected:
    // dir == -1 =>  FFT
    // dir ==  1 => IFFT
    BaseFft( int nPoints, int dir );

  private:
    int                   _nPoints;
    int                   _logPoints;
    double                _sqrtPoints;
    SArr<int>             _bitRev;

    SArr<Complex>         _x; // in-place fft array
    SArr< SArr<Complex> > _w; // exponentials
};

// Calculates forward FFT
class Fft : public BaseFft
{
  public:
    explicit Fft( int nPoints )
    : BaseFft( nPoints, -1 )
    {}
};

// Calculates inverse FFT
// Note: the result is NOT multiplied by normalizing factor 1/N
class Ifft : public BaseFft
{
  public:
    explicit Ifft( int nPoints )
    : BaseFft( nPoints, 1 )
    {}
};

} // namespace DSP

#endif // DSP_FFT_H
