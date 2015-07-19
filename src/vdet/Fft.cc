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

#include <cmath>
#include "Fft.h"

namespace DSP
{

BaseFft::BaseFft( int nPoints, int dir )
: _nPoints(nPoints),
  _logPoints(0),
  _sqrtPoints(::sqrt((double)_nPoints))
{
  // calculate binary log
  --nPoints;
  while ( nPoints != 0 )
  {
    nPoints >>= 1;
    ++_logPoints;
  }

  {
    SArr<int> temp( new int [_nPoints] );
    _bitRev = temp;
  }
  {
    SArr<Complex> temp( new Complex[_nPoints] );
    _x = temp;
  }
  {
    SArr< SArr<Complex> > temp( new SArr<Complex>[_logPoints+1] );
    _w = temp;
  }
  
  // Precompute complex exponentials
#ifndef PI
  double const PI = 2.0 * ::asin(1.0);
#endif
  int _2_l = 2;
  for ( int l = 1; l <= _logPoints; ++l )
  {
    {
      SArr<Complex> temp( new Complex[_nPoints] );
      _w[l] = temp;
    }
    
    for ( int i = 0; i < _nPoints; ++i )
    {
      double re =       ::cos( 2.0 * PI * i / _2_l );
      double im = dir * ::sin( 2.0 * PI * i / _2_l );
      _w[l][i] = Complex(re, im);
    }

    _2_l <<= 1;
  }
  
  // set up bit reverse mapping
  int rev = 0;
  int halfPoints = _nPoints / 2;
  for ( int i = 0; i < _nPoints - 1; ++i )
  {
    _bitRev[i] = rev;
    int mask = halfPoints;
    // add 1 backwards
    while ( rev >= mask )
    {
      rev -= mask; // turn off this bit
      mask >>= 1;
    }
    
    rev += mask;
  }
  
  _bitRev[_nPoints - 1] = _nPoints - 1;
}

} // namespace DSP

