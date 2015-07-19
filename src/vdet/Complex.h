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

#ifndef DSP_COMPLEX_H
#define DSP_COMPLEX_H

#include <cmath>

namespace DSP
{

class Complex
{
  public:
    Complex( void )
    : _re(0.0),
      _im(0.0)
    {}
    
    Complex( double re )
    : _re(re), 
      _im(0.0)
    {}
    
    Complex( double re, double im )
    : _re(re), 
      _im(im)
    {}
    
    double Re( void ) const 
    { return _re; }

    double Im( void ) const 
    { return _im; }

    void operator+=( Complex const& c )
    {
      _re += c._re;
      _im += c._im;
    }

    void operator-=( Complex const& c )
    {
      _re -= c._re;
      _im -= c._im;
    }

    void operator*=( Complex const& c )
    {
      double reT = c._re * _re - c._im * _im;
      _im = c._re * _im + c._im * _re;
      _re = reT;
    }

    Complex operator-( void ) 
    { return Complex( -_re, -_im ); }

    double Mod2( void ) const
    { return _re * _re + _im * _im; }

    double Mod( void ) const 
    { return ::sqrt( Mod2() ); }

    double Arg( void ) const
    { return ::atan2( _im, _re ); }
    
  private:
    double _re;
    double _im;
};

} // namespace DSP

#endif // DSP_COMPLEX_H
