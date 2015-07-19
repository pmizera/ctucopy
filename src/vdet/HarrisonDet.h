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

#ifndef VOICE_HARRISONDET_H
#define VOICE_HARRISONDET_H

/*
**

  Implementace Harrisonovy metody detekce rec - pauza 
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Volitelne parametry:

    NINITSEGMENTS - pocet inicializacnich segmentu
      
        Zpracovanim techto segmentu se nastavi pocatecni hodnoty
      pro odhad stredni hodnoty a stredni kvadraticke odchylky
      energie signalu, ze kterych se dale urci prahova hodnota
      energie pro detekci rec - pauza.

        Predpoklada se, ze tyto segmenty nebudou obsahovat recovy
      signal.

        Pro tyto segmenty je vysledek detekce vzdy false.


    P - adaptatacni konstanta pro aktualizaci prahove energie


  Iterator:

    Pro iterator musi byt definovany tyto operatory:

      operator++            (prefix i postfix verze)
      operator*             (dereference)
      operator== , resp. != (test rovnosti, resp. nerovnosti)

      Hodnota vracena dereferencnim operatorem * musi byt primo
    prevoditelna na typ double.


  Priklad pouziti:

    int const NSAMPLES = 256;    // pocet vzorku v segmentu
    short buf[NSAMPLES];

    short* begin = &buf[0];          // begin iterator
    short* end = &buf[0] + NSAMPLES; // end iterator

    HarrisonDetector detector;
    for ( int samplesRead = ReadSamples( buf, bufSize ); // Nacte 256 vzorku
          samplesRead == bufSize;                     // Mame cely segment?
          samplesRead = ReadSamples( buf, bufSize ) ) // Nacte dalsi vzorky
    {
      bool result = detector.Process( begin, end );
    }
**
*/

#include <cmath>
#include "Energy.h"
#include "Mean.h"
#include "StdDev.h"

namespace Voice
{

class HarrisonDetector
{
  public: // Class
    class Options
    {
      public:
        Options( int nInitSegments = 40, double p = 0.8, double q = 0.99, double za = 2.0 )
        : NINITSEGMENTS(nInitSegments),
          P(p),
          Q(q),
          ZA(za)
        {}

      public: // Data members
        int    NINITSEGMENTS; // number of initial segments
        double P; // adaptation constant for background energy updating
        double Q;
        double ZA;
    };

  public: // Methods
    HarrisonDetector( Options const& opt )
    : _opt(opt),
      _nSegmentsProcessed(0),
      _e0(0.0),
      _ep(0.0),
      _eMean(0.0),
      _eMean2(0.0),
      _eVar(0.0),
      _result(false)
    {}

    template <class SI>
    bool Process( SI begin, SI end )
    { 
      DSP::Energy en(begin,end);
      if ( _nSegmentsProcessed < _opt.NINITSEGMENTS )
      { // Initialization phase
        _result = false;

        if ( _nSegmentsProcessed < 2 )
        { _en[_nSegmentsProcessed] = en.GetEnergy();
        }
        else if ( _nSegmentsProcessed == 2 )
        { _eMean = DSP::Mean( &_en[0], &_en[2] ).GetMean();
          _eMean2 = ::pow( _eMean, 2.0 );
          double dev = DSP::StdDeviation( &_en[0], &_en[2], _eMean ).GetDeviation();
          _eVar = ::pow( dev, 2.0 );
          _e0 = _eMean;
          _ep = _e0 + _opt.ZA * dev;
        }
        else
        { double e = ::fabs(en.GetEnergy() - _e0);
          _e0 = _opt.P * _e0 + (1 - _opt.P) * en.GetEnergy();
          _eMean = _opt.Q * _eMean + (1 - _opt.Q) * e;
          _eMean2 = _opt.Q * _eMean2 + (1 - _opt.Q) * ::pow( e, 2.0 );
          _eVar = _eMean2 - ::pow( _eMean, 2.0 );
          _ep = _eMean + _opt.ZA * ::sqrt( _eVar );
        }

        ++_nSegmentsProcessed;
      }
      else
      { // Detection phase
        double e = ::fabs(en.GetEnergy() - _e0);
        _result = (e > _ep);

        if ( !_result )
        { // Update parameters
          _e0 = _opt.P * _e0 + (1 - _opt.P) * en.GetEnergy();
          _eMean = _opt.Q * _eMean + (1 - _opt.Q) * e;
          _eMean2 = _opt.Q * _eMean2 + (1 - _opt.Q) * ::pow( e, 2.0 );
          _eVar = _eMean2 - ::pow( _eMean, 2.0 );
          _ep = _eMean + _opt.ZA * ::sqrt( _eVar );
        }
      }

      return _result;
    }

    bool GetResult( void ) const
    { return _result; }

  private: // Data members
    Options  _opt;
    
    int      _nSegmentsProcessed;

    double   _en[3];
    double   _e0;
    double   _ep;
    double   _eMean;
    double   _eMean2;
    double   _eVar;
    
    bool     _result;
};

} // namespace Voice

#endif // VOICE_HARRISONDET_H
