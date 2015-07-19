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

#ifndef DSP_MEAN_H
#define DSP_MEAN_H

#include <cassert>

namespace DSP
{

class Mean
{
  public:
    Mean( void )
    : _mean(0.0)
    {}

    template <class SI>
    Mean( SI begin, SI end )
    : _mean(0.0)
    { Compute(begin,end); }

    template <class SI>
    double Compute( SI begin, SI end )
    {
      int n = 0;
      _mean = 0.0;
      for ( SI it = begin; it != end; ++it, ++n )
      { _mean += double(*it);
      }
      
      assert( n > 0 );
      _mean /= n;
      return _mean;
    }

    double GetMean( void ) const
    { return _mean; }

  private:
    double _mean;
};

} // namespace DSP

#endif // DSP_MEAN_H
