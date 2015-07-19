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

#ifndef SMARTPOINTER_H
#define SMARTPOINTER_H

#include <string.h>
#include <stdlib.h>
 
//
// Class SPtr (SmartPtr)
// ~~~~~~~~~~~~~~~~~~~~~~
template <class T>
class SPtr
{
  public:
	  // ~~ Contruction / Destruction~~
    SPtr( void )
    : _p(0) 
    {}
    
    explicit SPtr( T* p )
    : _p(p)
    {}

    ~SPtr( void )
    { delete _p; }

		// ~~ Contstruction from derived classes ~~~
		template <class U>
		SPtr( SPtr<U>& ptr )
		{ _p = ptr.release();	}
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Value semantics
    SPtr( SPtr<T>& ptr )
    { _p = ptr.release(); }


    SPtr<T>& operator=( SPtr<T>& ptr )
    {
      if ( _p != ptr._p )
      {
        delete _p;
        _p = ptr.release();
      }
			return *this;
    }

		// ~~ Assignment from derived classes ~~~
		template <class U>
		void operator=( SPtr<U>& ptr )
		{
		  if ( _p != ptr._p )
      {
        delete _p;
        _p = ptr.release();
      }
		}
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		template <class U>
		void down_cast( SPtr<U>& ptr )
		{
		  _p = dynamic_cast<T*>(ptr._p);
			if ( _p == 0 ) throw "Dynamic cast failed !";
			ptr.release();
		}

    T* release( void )
    {
      T* pTmp = _p;
      _p = 0;
      return pTmp;
    }

    T* operator->( void ) { return _p; }
    T const* operator->( void ) const { return _p; }
    T* get( void ) { return _p; }
    T const* get( void ) const { return _p; }
    
    T& operator* ( void ){ return *_p; }
    T const& operator*( void ) const { return *_p; }

		operator T*( void ) { return _p; }
    operator T const*( void ) const { return _p; }

  private:
    T* _p;
};

//
// Class SArr (SmartArray)
// ~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
class SArr
{
  public:
    SArr( void )
    : _p(0) 
    {}
    
    explicit SArr( T* p )
    : _p(p)
    {}

    ~SArr( void )
    { delete[] _p; }

    // Value semantics
    SArr( SArr<T>& ptr )
    { _p = ptr.release(); }


    void operator=( SArr<T>& ptr )
    {
      if ( _p != ptr._p )
      {
        delete[] _p;
        _p = ptr.release();
      }
    }

    T* release( void )
    {
      T* pTmp = _p;
      _p = 0;
      return pTmp;
    }

    T* get( void ) { return _p; }
    T const* get( void ) const { return _p; }
    
    T& operator* ( void ){ return _p[0]; }
    T const& operator*( void ) const { return _p[0]; }

    operator T*( void ) { return _p; }
    operator T const*( void ) const { return _p; }

    void qsort( int (*comp)( T const& p1, T const& p2 ), int from, int to )
    {
		  if ( from >= to ) return;
      int i = ::rand() % (to - from + 1) + from; // a random integer between from and to
      swap( from, i );
      int mid = from;
      for ( i = from + 1; i <= to; ++i )
      {  
			  if ( (*comp)( _p[i], _p[from] ) < 0 )
        {  
				  mid++;
          swap( mid, i );
        }
      }
      swap( from, mid );
      
     	qsort( comp, from, mid-1 );
      qsort( comp, mid+1, to );
    }

  private:
	  void swap( int i, int j )
		{
		  T temp = _p[i];
			_p[i] = _p[j];
			_p[j] = temp;
		}

  protected:
    T* _p;
};

//
// Class SStr (SmartString)
// ~~~~~~~~~~~~~~~~~~~~~~~~
class SStr : public SArr<char>
{
  public:
    SStr( void )
    {}

    explicit SStr( char const str[] ) // Makes copy of str
    : SArr<char>( new char[::strlen(str)+1] )
    { ::strcpy( _p, str ); }

    explicit SStr( int len ) // Excude terminating null
    : SArr<char>( new char[len+1] )
    { _p[0] = '\0'; }

    void Copy( char const str[] )
    {
      delete[] _p;
      _p = new char[::strlen(str)+1];
      ::strcpy( _p, str );
    }
};


#endif // SMARTPOINTER_H
