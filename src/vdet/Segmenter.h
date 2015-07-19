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

#ifndef DSP_SEGMENTER_H
#define DSP_SEGMENTER_H

#include <cassert>
#include "SmartPtr.h"

namespace DSP
{

/*

Trida _Input musi definovat metodu:

  int Read( void* buf, int bufSize );

ktera vraci pocet bytu nactenych do buf

*/

template <class _E, class _Input>
class Segmenter
{
  public:
    class Iterator
    {
      public: // Methods
        _E operator*( void ) const
        { return *_sample;
        }

        _E* operator++( void ) // Prefix version
        { ++_sample;
          if ( _sample == _end )
          { _wrapped = true;
            _sample = _begin;
          }
          return _sample;
        }

        _E* operator++( int ) // Postfix version
        { _E* temp = _sample;
          ++_sample;
          if ( _sample == _end )
          { _wrapped = true;
            _sample = _begin;
          }
          return temp;
        }

        bool operator==( Iterator const& b ) const
        { return (_wrapped || b._wrapped) && _sample == b._sample;
        }

        bool operator!=( Iterator const& b ) const
        { return !operator==(b);
        }

      private: // Methods
        Iterator( _E* begin, _E* end, _E* sample )
        : _begin(begin),
          _end(end),
          _sample(sample),
          _wrapped(false)
        {}

      private: // Data members
        _E*  _begin;
        _E*  _end;
        _E*  _sample;
        bool _wrapped;

      friend class Segmenter<_E, _Input>;
    };

  public:
    Segmenter( int nPoints, int step, _Input& input )
    : _nPoints(nPoints),
      _step(step),
      _input(input),
      _buf( new _E[_nPoints] ),
      _begin(&_buf[0]),
      _atEnd(false)
    { 
      assert( _step > 0 );
      assert( _step <= _nPoints );
      _atEnd = (_input.Read( _buf, _nPoints * sizeof(_E) ) != _nPoints * sizeof(_E));
    }

    bool AtEnd( void ) const
    { return _atEnd;
    }

    void Advance( void )
    { 
      if ( _begin + _step < &_buf[0] + _nPoints )
      { _atEnd = (_input.Read(_begin, _step * sizeof(_E)) != _step * sizeof(_E));
        _begin += _step;
      }
      else
      { int remainingEntries = (&_buf[0] + _nPoints) - _begin;
        _atEnd = (_input.Read(_begin, remainingEntries * sizeof(_E)) != remainingEntries * sizeof(_E));
        if ( !_atEnd && _step - remainingEntries > 0 )
        { _atEnd = (_input.Read(_buf, (_step - remainingEntries) * sizeof(_E)) != (_step - remainingEntries) * sizeof(_E));
        }
        _begin = &_buf[0] + _step - remainingEntries;
      }
    }

    Iterator Begin( void )
    { return Iterator( &_buf[0], &_buf[0] + _nPoints, _begin );
    }

    Iterator End( void )
    { return Iterator( &_buf[0], &_buf[0] + _nPoints, _begin );
      // Yes, the same as Begin, see equality operator for Iterator
      // (the _wrapped flag)
    }

  private:
    int const _nPoints;
    int const _step;
    _Input&   _input;
    SArr<_E>  _buf;
    _E*       _begin;
    bool      _atEnd;
};

} // namespace DSP

#endif // DSP_SEGMENTER_H
