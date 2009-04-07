/// \file ijk.txx
/// ijk templates defining general ijk objects, i.e. classes BOX and ERROR
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJK_
#define _IJK_

#include <iostream>
#include <sstream>
#include <vector>

namespace IJK {

  // **************************************************
  // TEMPLATE CLASS BOX
  // **************************************************

  template <class COORD_TYPE> class BOX {
  protected:
    int dimension;
    COORD_TYPE * min_coord;
    COORD_TYPE * max_coord;

    void Init();
    void FreeAll();
    
  public:
    BOX() { Init(); };
    BOX(const int dimension);
    ~BOX() { FreeAll(); };
    BOX(const BOX & box);                 // copy constructor
    const BOX & operator = (const BOX &); // copy assignment

    // get functions
    int Dimension() const { return(dimension); };
    COORD_TYPE MinCoord(const int d) const
    { return(min_coord[d]); };
    COORD_TYPE MaxCoord(const int d) const
    { return(max_coord[d]); };

    // set functions
    void SetDimension(const int d);
    void SetMinCoord(const int d, const COORD_TYPE c)
    { min_coord[d] = c; };
    void SetMaxCoord(const int d, const COORD_TYPE c)
    { max_coord[d] = c; };
    void SetAllMinCoord(const COORD_TYPE c);  // set all min coord to c
    void SetAllMaxCoord(const COORD_TYPE c);  // set all max coord to c
  };

  // **************************************************
  // TEMPLATE CLASS ARRAY
  // **************************************************

  /// Simple array class for creating static arrays
  /// Making this a class guarantees that when a program leaves
  /// the function where this class is declared, the array memory is freed.
  template <class ETYPE> class ARRAY {
  protected:
    ETYPE * element;

  public:
    template <class STYPE>
    ARRAY(const STYPE array_size) // constructor
    { element = new ETYPE[array_size]; };
    ~ARRAY();

    // get functions
    ETYPE * Ptr() { return(element); };
    const ETYPE * PtrConst() const { return(element); }
    template <class ITYPE>
    ETYPE & operator [] (const ITYPE i) { return(element[i]); };
    template <class ITYPE>
    ETYPE operator [] (const ITYPE i) const { return(element[i]); };

    // Free function
    void Free();
  };

  // **************************************************
  // ERROR CLASSES
  // **************************************************

  template < class T1, class T2, class T3 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3;
    return(os0.str());
  }

  template < class T1, class T2, class T3, class T4, class T5 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3, T4 a4, T5 a5)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5;
    return(os0.str());
  }

  template < class T1, class T2, class T3, class T4, class T5, 
	     class T6, class T7 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7;
    return(os0.str());
  }

  // error class
  class ERROR {

  protected:
    std::vector<std::string> msg;   // error messages

  public:
    ERROR(){};
    ERROR(const char * error_msg) { AddMessage(error_msg); };
    ERROR(const std::string & error_msg) { AddMessage(error_msg); };

    ERROR(const ERROR & error)                      // copy constructor
    {
      for (int i = 0; i < error.NumMessages(); i++)
	{ msg.push_back(error.Message(i)); };
    };

    const ERROR & operator = (const ERROR & right)  // copy assignment
    {
      if (&right != this) {
	for (int i = 0; i < right.NumMessages(); i++)
	  { msg.push_back(right.Message(i)); };
      }

      return(*this);
    };

    // get functions
    std::string Message(const int i) const
    {
      if (i >= 0 && i < NumMessages()) { return(msg[i]); }
      else { return(""); }
    }
    int NumMessages() const { return(msg.size()); };

    // get procedure error messages
    std::string ProcMessage(const std::string & procedure_name)
    {
      std::string error_msg = 
	"Error detected in procedure: " +  procedure_name + ".";
      return(error_msg);
    }
    std::string ProcMessage(const char * procedure_name)
    { return(ProcMessage(std::string(procedure_name))); }

    // set functions
    void AddMessage(const std::string & error_msg)
    { msg.push_back(error_msg); };
    void AddMessage(const char * error_msg)     
    { AddMessage(std::string(error_msg)); };
    void AddMessage(const char * part1, const int val, const char * part2)
    { AddMessage(compose_string(part1, val, part2)); };
    void AddMessage(const std::string & part1, const int val, 
		    const std::string & part2)
    { AddMessage(compose_string(part1, val, part2)); };
      std::ostringstream error_msg;
    void AddMessage(const std::string & part1, const int val, 
		    const char * part2)
    { AddMessage(compose_string(part1, val, part2)); };
    void AddMessage(const char * part1, const int val1, const char * part2,
		    const int val2, const char * part3)
    { AddMessage(compose_string(part1, val1, part2, val2, part3)); };
    void AddMessage(const char * part1, const int val1, const char * part2,
		    const int val2, const char * part3, const int val3,
		    const char * part4)
    { AddMessage(compose_string(part1, val1, part2, val2, part3, 
				val3, part4)); };
    void AddProcMessage(const char * procname)
    { AddMessage(ProcMessage(procname)); };
    void AddProcMessage(const std::string & procname)
    { AddMessage(ProcMessage(procname)); };
    void SetMessage(const int i, const std::string & error_msg)
    {
      if (i >= 0 && i < NumMessages())
	msg[i] = error_msg;
    };
    void SetMessage(const int i, const char * error_msg)
    { SetMessage(i, std::string(error_msg)); };
    void SetProcMessage(const int i, const char * procname)
    { SetMessage(i, ProcMessage(procname)); };
    void SetProcMessage(const int i, const std::string & procname)
    { SetMessage(i, ProcMessage(procname)); };
    void ClearAll() { msg.clear(); };

    ERROR & operator ()(const std::string & msg1)
      {
	AddMessage(msg1);
	return(*this);
      };

    ERROR & operator ()(const std::string & msg1, const std::string & msg2)
      {
	AddMessage(msg1);
	AddMessage(msg2);
	return(*this);
      }

    // print error messages to ostream out
    void Print(std::ostream & out)
    {
      for (int i = 0; i < NumMessages(); i++) {
	out << Message(i) << std::endl;
      }
    }

  };

  class PROCEDURE_ERROR:public ERROR {

  protected:
    std::string procedure_name;

  public:
    PROCEDURE_ERROR(const char * procedure_name)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
    };
    PROCEDURE_ERROR(const char * procedure_name, const char * error_msg)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
      AddMessage(error_msg);
    };
    PROCEDURE_ERROR(const std::string & procedure_name, 
		    const std::string & error_msg)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
      AddMessage(error_msg);
    };

  };

  // **************************************************
  // TEMPLATE CLASS BOX MEMBER FUNCTIONS
  // **************************************************

  template <class T> void BOX<T>::Init()
  {
    dimension = 0;
    min_coord = 0;
    max_coord = 0;
  }

  template <class T> void BOX<T>::FreeAll()
  {
    dimension = 0;
    if (min_coord != NULL) { delete [] min_coord; }
    if (max_coord != NULL) { delete [] max_coord; }
    min_coord = NULL;
    max_coord = NULL;
  }

  template <class T> BOX<T>::BOX(const int dimension)
  { Init(); SetDimension(dimension); };

  template <class COORD_TYPE> 
  void BOX<COORD_TYPE>::SetDimension(const int d)
  {
    FreeAll();
    Init();
    if (d <= 0) return;
    min_coord = new COORD_TYPE[d];
    max_coord = new COORD_TYPE[d];
    dimension = d;
  }


  template <class COORD_TYPE>
  void BOX<COORD_TYPE>::SetAllMinCoord(const COORD_TYPE c)
  {
    for (int d = 0; d < dimension; d++)
      SetMinCoord(d, c);
  }

  template <class COORD_TYPE>
  void BOX<COORD_TYPE>::SetAllMaxCoord(const COORD_TYPE c)
  {
    for (int d = 0; d < dimension; d++)
      SetMaxCoord(d, c);
  }

  // copy constructor template
  template <class COORD_TYPE> 
  BOX<COORD_TYPE>::BOX(const BOX<COORD_TYPE> & box)
  {
    Init();
    SetDimension(box.Dimension());
    for (int d = 0; d < box.Dimension(); d++) {
      SetMinCoord(d, box.MinCoord(d));
      SetMaxCoord(d, box.MaxCoord(d));
    }
  }

  // copy assigment template
  template <class COORD_TYPE> 
  const BOX<COORD_TYPE> & BOX<COORD_TYPE>::operator = 
    (const BOX<COORD_TYPE> & right)
  {
    if (&right != this) {         // avoid self-assignment
      FreeAll();
      SetDimension(right.Dimension());
      for (int d = 0; d < right.Dimension(); d++) {
	SetMinCoord(d, right.MinCoord(d));
	SetMaxCoord(d, right.MaxCoord(d));
      }
    }
    return *this;
  }

  // **************************************************
  // TEMPLATE CLASS ARRAY MEMBER FUNCTIONS
  // **************************************************

  template <class ETYPE> ARRAY<ETYPE>::~ARRAY()
  {
    Free();
  }

  template <class ETYPE> void ARRAY<ETYPE>::Free()
  {
    delete [] element;
    element = NULL;
  }

}

#endif
