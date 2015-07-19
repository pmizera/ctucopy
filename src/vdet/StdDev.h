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

#ifndef DSP_STDDEV_H
#define DSP_STDDEV_H

#include <cassert>
#include <cmath>
#include "Mean.h"

namespace DSP
{

class StdDeviation
{
  public:
    StdDeviation( void )
    : _dev(0.0)
    {}

    template <class SI>
    StdDeviation( SI begin, SI end )
    : _dev(0.0)
    { Compute(begin,end); }

    template <class SI>
    StdDeviation( SI begin, SI end, double mean )
    : _dev(0.0)
    { Compute(begin,end,mean); }


    template <class SI>
    double Compute( SI begin, SI end )
    { return Compute( begin, end, Mean( begin, end ).GetMean() );
    }

    template <class SI>
    double Compute( SI begin, SI end, double mean )
    {
      int n = 0;
      _dev = 0.0;
      for ( SI it = begin; it != end; ++it, ++n )
      { _dev += ::pow( ::fabs( double(*it) - mean ), 2.0 );
      }

      assert( n > 1 );
      _dev = ::sqrt( _dev / (n - 1) );
      return _dev;
    }

    double GetDeviation( void ) const
    { return _dev; }

  private:
    double _dev;
};

} // namespace DSP

#endif // DSP_STDDEV_H
