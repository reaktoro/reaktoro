//-------------------------------------------------------------------
// $Id: gdatastream.cpp 771 2012-12-13 13:07:43Z kulik $
//
/// \file gdatastream.cpp
/// Implementation of stream binary file operations extended for endianness
/// (e.g. for compatibility between Intel- and old Mac processors)
//
// Copyright (c) 1996-2012 A.Rysin, S.Dmytriyeva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <stdint.h>
using namespace std;

#ifdef __linux__
#include <endian.h>
#elif defined(__APPLE__)
#include <machine/endian.h>
#define __BYTE_ORDER BYTE_ORDER
#define __BIG_ENDIAN BIG_ENDIAN
#else
//#warning "other plaftorm - considering little endian"
#define __BIG_ENDIAN 4321
#define __BYTE_ORDER 1234
#endif

#include "gdatastream.h"

//--------------------------------------------------------------------

#ifndef __BYTE_ORDER
# error "Error: __BYTE_ORDER not defined\n";
#endif

inline short SWAP(short x) { 
    return (((x>>8) & 0x00ff) | ((x<<8) & 0xff00)); 
}

inline int32_t SWAP(int32_t x) {
    int32_t x_new;
    char* xc_new = (char*)&x_new;
    char* xc = (char*)&x;
    xc_new[0] = xc[3];
    xc_new[1] = xc[2];
    xc_new[2] = xc[1];
    xc_new[3] = xc[0];
    return x_new;
}

/*inline int SWAP(int x) {
    int x_new;
    char* xc_new = (char*)&x_new;
    char* xc = (char*)&x;
    xc_new[0] = xc[3];
    xc_new[1] = xc[2];
    xc_new[2] = xc[1];
    xc_new[3] = xc[0];
    return x_new;
}*/
/*inline int SWAP(long x) {
    long x_new;
    char* xc_new = (char*)&x_new;
    char* xc = (char*)&x;
    xc_new[0] = xc[3];
    xc_new[1] = xc[2];
    xc_new[2] = xc[1];
    xc_new[3] = xc[0];
    return x_new;
}*/


/*
inline long SWAP(long x) {
    long x_new;
    char* xc_new = (char*)&x_new;
    char* xc = (char*)&x;
    xc_new[0] = xc[7];
    xc_new[1] = xc[6];
    xc_new[2] = xc[5];
    xc_new[3] = xc[4];
    xc_new[4] = xc[3];
    xc_new[5] = xc[2];
    xc_new[6] = xc[1];
    xc_new[7] = xc[0];
    return x_new;
}
*/
inline float SWAP(float x) {
    float x_new;
    char* xc_new = (char*)&x_new;
    char* xc = (char*)&x;
    xc_new[0] = xc[3];
    xc_new[1] = xc[2];
    xc_new[2] = xc[1];
    xc_new[3] = xc[0];
    return x_new;
}

inline double SWAP(double x) {
    double x_new;
    char* xc_new = (char*)&x_new;
    char* xc = (char*)&x;
    xc_new[0] = xc[7];
    xc_new[1] = xc[6];
    xc_new[2] = xc[5];
    xc_new[3] = xc[4];
    xc_new[4] = xc[3];
    xc_new[5] = xc[2];
    xc_new[6] = xc[1];
    xc_new[7] = xc[0];
    return x_new;
}


GemDataStream::GemDataStream( gstring& aPath, ios::openmode aMod  ):
        mod( aMod ),
        Path( aPath ),
    //    byteorder( LittleEndian ),
        ff(aPath.c_str(), aMod)
{
    setByteOrder(LittleEndian);
    ErrorIf( !ff.good(), Path.c_str(), "Fileopen error");
}

GemDataStream::~GemDataStream()
{
}


void GemDataStream::setByteOrder( int bo )
{
    byteorder = bo;

#if __BYTE_ORDER == __BIG_ENDIAN
	swap = (byteorder == LittleEndian);
#warning "Compiling for BIG ENDIAN architecture!"
#else
	swap = (byteorder == BigEndian);
#endif
//    cerr << "GemDataStream::swap == " << swap << endl;
}

// NOTE: these functions better to write as a templates!!

GemDataStream &GemDataStream::operator>>( char &i )
{
    ff.read((char*)&i, sizeof(char));
    return *this;
}

GemDataStream &GemDataStream::operator>>( short &i )
{
    ff.read((char*)&i, sizeof(short));
    if( swap ) i = SWAP(i);
    return *this;
}

GemDataStream &GemDataStream::operator>>( int &i_ )
{
    //ff.read((char*)&i, sizeof(int));
    int32_t i=(int32_t)i_;
    ff.read((char*)&i, sizeof(int32_t));
    if( swap ) i = SWAP(i);
    i_ = (int)i;
    return *this;
}

GemDataStream &GemDataStream::operator>>( long &i_ )
{
    //ff.read((char*)&i, sizeof(long));
    int32_t i=(int32_t)i_;
    ff.read((char*)&i, sizeof(int32_t));
    if( swap ) i = SWAP(i);
    i_ = (long)i;
    return *this;
}

GemDataStream &GemDataStream::operator>>( float &f )
{
    ff.read((char*)&f, sizeof(float));
    if( swap ) f = SWAP(f);
    return *this;
}

GemDataStream &GemDataStream::operator>>( double &f )
{
    ff.read((char*)&f, sizeof(double));
    if( swap ) f = SWAP(f);
    return *this;
}


// NOTE: these functions are inefficient !!!
// it's faster to read the whole array
// and then loop through it to reverse byte order!!!


GemDataStream &GemDataStream::operator<<( char i )
{
    ff.write(&i, sizeof(char));
    return *this;
}

GemDataStream &GemDataStream::operator<<( short i )
{
    if( swap ) i = SWAP(i);
    ff.write((char*)&i, sizeof(short));
    return *this;
}

GemDataStream &GemDataStream::operator<<( int i_ )
{
    int32_t i=(int32_t)i_;
    if( swap ) i = SWAP(i);
    ff.write((char*)&i, sizeof(int32_t));
    return *this;
}

GemDataStream &GemDataStream::operator<<( long i_ )
{
    int32_t i=(int32_t)i_;
    if( swap ) i = SWAP(i);
    ff.write((char*)&i, sizeof(int32_t));
    return *this;
}


GemDataStream &GemDataStream::operator<<( float f )
{
    if( swap ) f = SWAP(f);
    ff.write((char*)&f, sizeof(float));
    return *this;
}


GemDataStream &GemDataStream::operator<<( double f )
{
    if( swap ) f = SWAP(f);
    ff.write((char*)&f, sizeof(double));
    return *this;
}

void GemDataStream::readArray( char* arr, int size )
{
  if( !arr )
    return;
    
  ff.read(arr, size);
//  for(int ii=0; ii<size; ii++)
//   *this >> arr[ii];
}

void GemDataStream::readArray( short* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this >> arr[ii];
}

void GemDataStream::readArray( int* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this >> arr[ii];
}

void GemDataStream::readArray( long* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this >> arr[ii];
}

void GemDataStream::readArray( float* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this >> arr[ii];
}

void GemDataStream::readArray( double* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this >> arr[ii];
}

void GemDataStream::writeArray( char* arr, int size )
{
  if( !arr )
    return;

  ff.write(arr, size);
//  for(int ii=0; ii<size; ii++)
//   *this << arr[ii];
}

void GemDataStream::writeArray( short* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this << arr[ii];
}

void GemDataStream::writeArray( int* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this << arr[ii];
}

void GemDataStream::writeArray( long* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this << arr[ii];
}

void GemDataStream::writeArray( float* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this << arr[ii];
}
void GemDataStream::writeArray( double* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this << arr[ii];
}
// gdatastream.cpp

