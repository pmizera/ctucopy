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

#ifndef DSP_ENERGY_H
#define DSP_ENERGY_H

#include <cmath>

namespace DSP
{

class Energy
{
  public:
    Energy( void )
    : _en(0.0)
    {}

    template <class SI>
    Energy( SI begin, SI end )
    : _en(0.0)
    { Compute(begin,end); }

    template <class SI>
    double Compute( SI begin, SI end )
    { _en = 0.0;
      for ( SI it = begin; it != end; ++it )
      { _en += ::pow( double(*it), 2.0 );
      }
      return _en;
    }

    double GetEnergy( void ) const
    { return _en; }

  private:
    double _en;
};

} // namespace DSP

#endif // DSP_ENERGY_H
