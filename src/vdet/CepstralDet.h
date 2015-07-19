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

#ifndef VOICE_CEPSTRALDET_H
#define VOICE_CEPSTRALDET_H

#include <cassert>
#include <cmath>
#include "SmartPtr.h"
#include "Burg.h"
#include "Fft.h"

/*

class CepstrumEstimator
{
  public:
    CepstrumEstimator( int nPoints, int nCoefs );

    template <class SI>
    void Process( SI begin, SI end ); // computes the first !! nCoefs !! 
                                      // cepstral coefficients from
                                      // !! nPoints !! samples

    double operator[]( int i ) const; // returns i-th cepstral coefficient

    iterator Begin();                 // for iteration through
    iterator End();                   // the !! nCoefs !! coefficients
                                      // !! NOT nPoints !!
};

*/

namespace Voice
{

class CepstralDistance
{
  public:
    CepstralDistance( void )
    : _d(0.0)
    {}

    template <class SI1, class SI2>
    CepstralDistance( SI1 begin1, SI1 end1, SI2 begin2, SI2 end2 )
    : _d(0.0)
    { Compute( begin1, end1, begin2, end2 );
    }

    template <class SI1, class SI2>
    double Compute( SI1 begin1, SI1 end1, SI2 begin2, SI2 end2 )
    {
      double sum = 0.0;
      SI1 it1 = begin1;
      SI2 it2 = begin2;
      for ( ++it1, ++it2; // Skipping the first coefs
            it1 != end1 && it2 != end2;
            ++it1, ++it2 )
      {
        sum += (double(*it1) - double(*it2)) * (double(*it1) - double(*it2));
      }

      assert( it1 == end1 && it2 == end2 ); // The vectors aren't 
                                            // of same length

      _d = 4.3429 * ::sqrt(2 * sum); 
      return _d;
    }

    double GetDistance( void ) const
    { return _d;
    }

  private:
    double _d;
};

template <class _CepstrumEstimator>
class CepstralDetector
{
  public: // Class
    class Options
    {
      public:
        Options( int nInitSegments = 40, int nCoefs = 10, double p = 0.8, double q = 0.97 )
        : NINITSEGMENTS(nInitSegments),
          NCOEFS(nCoefs),
          P(p),
          Q(q)
        {}

      public:
        int    NINITSEGMENTS; // how many segments are guaranteed to be voiceless
        int    NCOEFS;        // number of cepstral coeffs for computing cepstral distance
        double P; // adaptation constant for cepstr. coefs
        double Q; // adaptation constant for threshold params
    };

  public: // Methods
    CepstralDetector( int nPoints, Options const& opt = Options() )
    : _opt(opt),
      _nPoints(nPoints),
      _han( new double [_nPoints] ),
      _x( new double[_nPoints] ),
      _c0( new double[_opt.NCOEFS] ),
      _ci( _nPoints, _opt.NCOEFS ),
      _dMean(0.0),
      _dMean2(0.0),
      _dVar(0.0),
      _threshold(0.0),
      _result(false),
      _nSegmentsProcessed(0)
    { // Compute Hann window coefficients
      double m = 2 * 3.141592653 / _nPoints;
      for ( int i = 0; i < _nPoints; ++i )
      { _han[i]= 0.5 * (1 - ::cos(m * i));
      }
    }

    template <class SI>
    bool Process( SI begin, SI end )
    {
      _result = false;
      // Window the samples
      int i = 0;
      SI it = begin;
      for ( ; it != end && i < _nPoints; ++it, ++i )
      { _x[i] = _han[i] * double(*it);
      }

      assert( it == end && i == _nPoints ); // The input vector doesn't have
                                            // _nPoints samples

      _ci.Process( &_x[0], &_x[0] + _nPoints );

      if ( _nSegmentsProcessed == 0 )
      {
        for ( int i = 0; i < _opt.NCOEFS; ++i )
        { _c0[i] = _ci[i];
        }
      } 
      else if ( _nSegmentsProcessed == 1 ) 
      {
	      for ( int i = 0; i < _opt.NCOEFS; ++i )
        { _c0[i] = (_c0[i] + _ci[i]) / 2.0;
        }

        double dist = CepstralDistance( _ci.Begin(), 
                                        _ci.End(), 
                                        &_c0[0], 
                                        &_c0[0] + _opt.NCOEFS ).GetDistance();

        _dMean = dist;
        _dMean2 = dist * dist;
        _threshold = _dMean;// + 2.0 * ::sqrt(_dVar); // _dVar == 0
      }
      else 
      {
        double dist = CepstralDistance( _ci.Begin(), 
                                        _ci.End(), 
                                        &_c0[0], 
                                        &_c0[0] + _opt.NCOEFS ).GetDistance();

        _result = (_nSegmentsProcessed > _opt.NINITSEGMENTS && dist >= _threshold);
        if ( !_result )
        { // no voice
          for ( int i = 0; i < _opt.NCOEFS; ++i )
          { _c0[i] = _opt.P * _c0[i] + (1 - _opt.P) * _ci[i];
          }

		      _dMean = _opt.Q * _dMean + (1 - _opt.Q) * dist;
		      _dMean2 = _opt.Q * _dMean2 + (1 - _opt.Q) * dist * dist;
		      _dVar = _dMean2 - _dMean * _dMean;
		      _threshold = _dMean + 2.0 * ::sqrt(_dVar);
        }
      }

      ++_nSegmentsProcessed; // !! may overflow => new initialization !!
      return _result;
    }

    bool GetResult( void ) const
    { return _result;
    }

  private:
    Options            _opt;
    int const          _nPoints;
    SArr<double>       _han;     // Hann window coefficients
    SArr<double>       _x;       // input buffer
    SArr<double>       _c0;

    _CepstrumEstimator _ci;

    double             _dMean;
    double             _dMean2;
    double             _dVar;
    double             _threshold;

    bool               _result;

    int                _nSegmentsProcessed;
};


class BurgCepstrumEstimator
{
  public: // Typedefs
    typedef DSP::Burg2Cepstrum::const_iterator const_iterator;

  public: // Methods
    BurgCepstrumEstimator( int nPoints, int nCoefs )
    : _a( nPoints, nCoefs ),
      _c( nCoefs )
    {}

    template <class SI>
    void Process( SI begin, SI end )
    { _a.Process(begin,end);
      _c.Process(_a);
    }

    double operator[]( int i )
    { return _c[i];
    }

    const_iterator Begin( void ) const
    { return _c.Begin();
    }

    const_iterator End( void ) const
    { return _c.End();
    }

  private:
    DSP::Burg          _a;
    DSP::Burg2Cepstrum _c;
};


class RealCepstrumEstimator
{
  public: // Typedefs
    typedef double const* const_iterator;

  public:
    // Note: nPoints must be power of 2 (because of FFT)
    RealCepstrumEstimator( int nPoints, int nCoefs )
    : _nPoints(nPoints),
      _nCoefs(nCoefs),
      _fft(_nPoints),
      _ifft(_nPoints),
      _c( new double[_nPoints] )
    {}

    double operator[]( int i ) const
    { 
      assert( i >= 0 );
      assert( i < _nCoefs );
      return _c[i];
    }

    const_iterator Begin( void ) const
    { return &_c[0]; 
    }

    const_iterator End( void ) const
    { return &_c[0] + _nCoefs; // !! NOT  + _nPoints !!
    }

    template<class SI>
    void Process( SI begin, SI end )
    {
      _fft.Process( begin, end );
      for ( int i = 0; i < _nPoints; ++i )
      { _c[i] = ::log( _fft.Mod(i) );
      }
      _ifft.Process( &_c[0], &_c[0] + _nPoints );
      for ( int i = 0; i < _nCoefs; ++i )
      { _c[i] = _ifft.Re(i) / _nPoints;
      }
    }

  private: // Data members
    int          _nPoints;
    int          _nCoefs;
    DSP::Fft     _fft;
    DSP::Ifft    _ifft;
    SArr<double> _c;
};

} // namespace Voice

#endif // VOICE_CEPSTRALDET_H
