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

#ifndef DSP_BURG_H
#define DSP_BURG_H

#include <cassert>
#include <cmath>
#include "SmartPtr.h"
#include "Energy.h"

namespace DSP
{

class Burg
{
  public: // Typedefs
    typedef double const* const_iterator;

  public: // Methods
    Burg( int nPoints, int nCoefs )
    : _nPoints(nPoints),
      _nCoefs(nCoefs),
      _ef( new double[_nPoints] ),
      _eb( new double[_nPoints] ),
      _efOld( new double[_nPoints] ),
      _ebOld( new double[_nPoints] ),
      _a( new double[_nCoefs] ),
      _aa( new double[_nCoefs] ),
      _rc( new double[_nCoefs] ),
      _alpha(0.0)
    { _a[0] = 1.0; // Always
    }

    template <class SI>
    void Process( SI begin, SI end )
    {
      int i; SI it;
      for ( i = 0, it = begin; it != end && i < _nPoints; ++it, ++i )
      { _ef[i] = _eb[i] = double(*it);
      }

      assert( it == end && i == _nPoints );

      // inic. vektoru dopredne a zpetne chyby predikce
      _alpha = Energy(begin,end).GetEnergy() / _nPoints;

      for ( int ik = 1; ik < _nCoefs; ++ik )
      { // vypocet parcialni korelace a koeficientu odrazu v sekci ik
        double num = 0.0;
        double den = 0.0;
        for ( int i = ik; i < _nPoints; ++i )
        { den += _ef[i] * _ef[i] + _eb[i - 1] * _eb[i - 1];
          num += _ef[i] * _eb[i - 1];
        }

        num *= 2.0;
        assert( den != 0.0 );
        _a[ik] = _rc[ik] = - num / den;
        _alpha *= 1 - _rc[ik] * _rc[ik];

        // filtrace chybovych signalu
        for ( int i = 0; i < _nPoints; ++i )
        { _efOld[i] = _ef[i];
          _ebOld[i] = _eb[i];
        }

        for ( int i = 1; i < _nPoints; ++i )
        { _ef[i] = _ef[i] + (_rc[ik] * _ebOld[i - 1]);
          _eb[i] = _ebOld[i - 1] + (_rc[ik] * _efOld[i]);
        }

        for( int i = 1; i < ik; ++i )
        { _a[i] = _aa[i] + _rc[ik] * _aa[ik - i];
        }

        for ( int i = 1; i <= ik; ++i )
        { _aa[i] = _a[i];
        }
      }
    }

    double GetAlpha( void ) const
    { return _alpha;
    }

    double operator[]( int i ) const
    { assert( i >= 0 );
      assert( i < _nCoefs );
      return _a[i];
    }

    const_iterator Begin( void ) const
    { return &_a[0];
    }

    const_iterator End( void ) const
    { return &_a[0] + _nCoefs;
    }

  private:
    int const    _nPoints;
    int const    _nCoefs;

    SArr<double> _ef;
    SArr<double> _eb;
    SArr<double> _efOld;
    SArr<double> _ebOld;

    SArr<double> _a;
    SArr<double> _aa;
    SArr<double> _rc;
    double       _alpha;
};

class Burg2Cepstrum
{
  public: // Typedefs
    typedef double const* const_iterator;

  public: // Methods
    explicit Burg2Cepstrum( int nCoefs )
    : _nCoefs(nCoefs),
      _c( new double[_nCoefs] )
    {}

    void Process( Burg const& burg )
    {
      for( int n = 1; n < _nCoefs; ++n )
      {
        double sum = 0.0;
        for ( int k = 1; k < n; ++k )
        { sum += (n-k) * _c[n-k] * burg[k];
        }
        _c[n]= - burg[n] - sum / n;
      }
      _c[0] = ::log( burg.GetAlpha() );
    }

    double operator[]( int i ) const
    { assert( i >= 0 );
      assert( i < _nCoefs );
      return _c[i];
    }

    const_iterator Begin( void ) const
    { return &_c[0];
    }

    const_iterator End( void ) const
    { return &_c[0] + _nCoefs;
    }

  private:
    int const    _nCoefs;
    SArr<double> _c;
};

} // namespace DSP

#endif // DSP_BURG_H
