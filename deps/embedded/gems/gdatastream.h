//-------------------------------------------------------------------
// $Id: gdatastream.h 771 2012-12-13 13:07:43Z kulik $
/// \file gdatastream.h
/// Contains definition of the GemDataStream class for processing
/// binary data streams on platforms with different endianness.
//
/// \class  GemDataStream  gdatastream.h
/// Stream binary file operations extended for endianness
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

#ifndef _gemdatastream_h_
#define _gemdatastream_h_

#include <fstream>
#include "verror.h"

class GemDataStream				// data stream class
{

    ios::openmode mod;
    gstring Path;

    int swap;
    int	byteorder;
    fstream ff;

public:
//    GemDataStream( fstream& ff  );
    GemDataStream( ) {    setByteOrder(LittleEndian); };
    GemDataStream( gstring& aPath, ios::openmode aMod  );
    virtual ~GemDataStream();

    const gstring& GetPath() const
    {
        return Path;
    }
    
//    bool	 atEnd() const;
//    bool	 eof() const;

    enum ByteOrder { BigEndian, LittleEndian };
    int	 byteOrder() const { return byteorder; };
    void setByteOrder( int );

    filebuf* rdbuf() { return ff.rdbuf(); }
    streamsize gcount() { return ff.gcount(); }
    istream& getline(char* s, streamsize n, char delim) { return ff.getline(s, n, delim); }
    void close() { ff.close(); }
    void put(char ch) { ff.put(ch); }
    istream& get(char& ch) { return ff.get(ch); }
    void sync() { ff.sync(); }
    bool good() { return ff.good(); }
    void clear() { ff.clear(); }
    void flush() { ff.flush(); }
    size_t tellg() { return ff.tellg(); }
    void open(const char* filename, ios::openmode mode) { ff.open(filename, mode); }
    ostream& seekp(size_t pos, ios_base::seekdir dir) { return ff.seekp(pos, dir); }
    istream& seekg(size_t pos, ios_base::seekdir dir) { return ff.seekg(pos, dir); }

    GemDataStream &operator>>( char &i );
    GemDataStream &operator>>( unsigned char &i ) { return operator>>((char&)i); }
    GemDataStream &operator>>( signed char &i ) { return operator>>((char&)i); }
    GemDataStream &operator>>( short &i );
    GemDataStream &operator>>( unsigned short &i ) { return operator>>((short&)i); }
    GemDataStream &operator>>( int &i );
    GemDataStream &operator>>( unsigned int &i ) { return operator>>((int&)i); }
    GemDataStream &operator>>( long &i );
    GemDataStream &operator>>( unsigned long &i ) { return operator>>((long&)i); }
    GemDataStream &operator>>( float &f );
    GemDataStream &operator>>( double &f );
//    GemDataStream &operator>>( char *&str );

    GemDataStream &operator<<( char i );
    GemDataStream &operator<<( unsigned char i ) { return operator<<((char) i); }
    GemDataStream &operator<<( signed char i ) { return operator<<((char) i); }
    GemDataStream &operator<<( short i );
    GemDataStream &operator<<( unsigned short i ) { return operator<<((short) i); }
    GemDataStream &operator<<( int i );
    GemDataStream &operator<<( unsigned int i ) { return operator<<((int) i); }
    GemDataStream &operator<<( long i );
    GemDataStream &operator<<( unsigned long i ) { return operator<<((long) i); }
    GemDataStream &operator<<( float f );
    GemDataStream &operator<<( double f );
//    GemDataStream &operator<<( const char *str );

    void readArray( char* arr, int size );
    void readArray( short* arr, int size );
    void readArray( int* arr, int size );
    void readArray( long* arr, int size );
    void readArray( float* arr, int size );
    void readArray( double* arr, int size );

    void writeArray( char* arr, int size );
    void writeArray( short* arr, int size );
    void writeArray( int* arr, int size );
    void writeArray( long* arr, int size );
    void writeArray( float* arr, int size );
    void writeArray( double* arr, int size );

#ifndef IPMGEMPLUGIN

    template <class T> void writeArray( T* arr, int size );
    template <class T> void readArray( T* arr, int size );
#endif
};

#ifndef IPMGEMPLUGIN

template <class T>
void GemDataStream::writeArray( T* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this << arr[ii];
}

template <class T>
void GemDataStream::readArray( T* arr, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
   *this >> arr[ii];
}

#endif
#endif
