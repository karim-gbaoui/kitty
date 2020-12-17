/* kitty: C++ truth table library
 * Copyright (C) 2017-2020  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file threshold_identification.hpp
  \brief Threshold logic function identification

  \author CS-472 2020 Fall students
*/

#pragma once

#include <vector>
#include <fstream>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"
#include "isop.hpp"
#include "kitty.hpp"
#include "operations.hpp"
#include "cube.hpp"
#include "properties.hpp"
#include "static_truth_table.hpp"
#include "dynamic_truth_table.hpp"
#include <sstream>
#include <string>


namespace kitty
{

/*! \brief Threshold logic function identification

  Given a truth table, this function determines whether it is a threshold logic function (TF)
  and finds a linear form if it is. A Boolean function is a TF if it can be expressed as

  f(x_1, ..., x_n) = \sum_{i=1}^n w_i x_i >= T

  where w_i are the weight values and T is the threshold value.
  The linear form of a TF is the vector [w_1, ..., w_n; T].

  \param tt The truth table
  \param plf Pointer to a vector that will hold a linear form of `tt` if it is a TF.
             The linear form has `tt.num_vars()` weight values and the threshold value
             in the end.
  \return `true` if `tt` is a TF; `false` if `tt` is a non-TF.
*/

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
  std::vector<int64_t> linear_form;
  auto numvars = tt.num_vars();
  TT tt0=tt;
  for ( auto i = 0u; i < numvars; i++ )
  {

    bool neg_unate = is_negative_unate_in_i(tt,i);
    if ( is_positive_unate_in_i(tt,i) + neg_unate == 0)
    {
     return false;
    }
    if (neg_unate){
	flip_inplace(tt0,i);
	}
    
  }

    dump_lp( tt0, "thresholdfunction.lp" );
    std::system( "lp_solve thresholdfunction.lp  > solution.txt" );	
    std::ifstream fin( "solution.txt", std::ifstream::in );
    

    /* parsing */
    std::string line, obj;
    std::getline( fin, line ); /* first line is empty */
    if (line == "This problem is infeasible"){
	return false;
    } 

    std::getline( fin, line ); /* first line is empty */
    std::getline( fin, obj ); /* second line is the value of objective function */
    std::getline( fin, line ); /* third line is empty */
    //std::getline( fin, line ); /* fourth line is useless */
    
    while ( std::getline( fin, line ) )
    {
	std::stringstream ss;
        ss<< line ;
	std::string temp_str;
	int temp_int;
	while(!ss.eof()) {
		ss >> temp_str;
	if(std::stringstream(temp_str) >> temp_int) { 
	
	 linear_form.push_back(uint64_t(temp_int));
         }
	temp_str = "";
	}
    }
    
   std::ofstream MyFile("filename.txt");
  
   for (auto i = 0 ; i<linear_form.size()-1; i++){
	if (is_negative_unate_in_i(tt,i)){
	linear_form[i] = -linear_form[i];
	linear_form[linear_form.size()-1] = linear_form[linear_form.size()-1] +linear_form[i];
		} 		
	}
   for (auto i = 0 ; i<linear_form.size(); i++){
	MyFile << linear_form[i] <<std::endl;
	}
   MyFile.close();

	
  /* if tt is TF: */
  /* push the weight and threshold values into `linear_form` */
  if ( plf )
  {
    *plf = linear_form;
  }
  return true;
  
}

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
void print_lp( const TT& tt, std::ostream& os = std::cout )
{
  /* the objective function */
  os << "min: ";
  for ( uint32_t i = 0; i < tt.num_vars(); i++ )
  {
    
    os << "+ "
       << "weight" << int(i);
  }
  os << " +Thresh;" << std::endl;

  /* onset constraints */
  std::vector<std::vector<uint8_t>> on_const = on_set_constraints( tt );
  for ( auto a : on_const )
  {
    if (a.empty()){ os<< "+ 0 " << " => " << " Thresh;"<<std::endl;
	}
    else{
    for ( auto b : a )
    {
      os << " +"
         << "weight" << int(b);
    }
    os << ">="
       << " Thresh;" << std::endl;
       }
    }

  /* offset constraints */
  std::vector<std::vector<uint8_t>> off_const = off_set_constraints( tt );
  for ( auto a : off_const )
  {
    if (a.empty()){ os<< "+ 0 " << " <= " << " Thresh-1;"<<std::endl;}
    else{
    for ( auto b : a )
    {
      os << "+ " <<"weight"<<int(b);
    }
    os << " <= "<< " Thresh-1;" << std::endl;
    }
  }

  /* variable type declaration */
  os << "int ";
  for ( uint8_t i = 0; i < tt.num_vars(); i++ )
  {
    os << "weight" << int(i) << ", ";
  }

  os << "Thresh"<< ";" << std::endl;

}

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
void dump_lp( const TT& tt, std::string const& filename = "thresholdfunction.lp" ) 
{
  std::ofstream fout(filename, std::ofstream::out );
  print_lp(tt, fout);
}


template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
std::vector<std::vector<uint8_t>> on_set_constraints( const TT& tt )
{

    std::vector<cube> on_set = isop(tt);
    std::vector<std::vector<uint8_t>> onset_const;
    for ( auto a : on_set )
    {
      std::vector<uint8_t> intermediate;
      for ( uint8_t i = 0; i < tt.num_vars(); i++ )
      {
        if ( a.get_mask(i)==1 && a.get_bit(i)==1 )
        {
          intermediate.push_back( i );
        }
      }
      onset_const.push_back( intermediate );
    }
    return onset_const;
}

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
TT operator~(const TT& tt)
{
    TT inter = tt;
    inter.num_vars = tt.num_vars;
    inter.bits = ~ tt.bits;	
    return inter;
}


template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
std::vector<std::vector<uint8_t>> off_set_constraints( const TT& tt )
{

  std::vector<cube> off_set = isop( operator~(tt));
  std::vector<std::vector<uint8_t>> offset_const;
  for ( auto a : off_set )
  {
    std::vector<uint8_t> intermediate;
    for ( uint8_t i = 0; i < tt.num_vars(); i++ )
    {
      if ( !a.get_mask(i) || (a.get_mask(i) && a.get_bit(i)) )
      {
        intermediate.push_back( i );
      }
    }
    offset_const.push_back( intermediate );
  }
  return offset_const;
}

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_negative_unate_in_i( const TT& tt, uint8_t i )
{
  auto const tt0 = cofactor0( tt, i );
  auto const tt1 = cofactor1( tt, i );
  for ( auto bit = 0; bit < ( 2 << ( tt.num_vars() - 1 ) ); bit++ )
  {
    if ( get_bit( tt0, bit ) >= get_bit( tt1, bit ) )
    {
      continue;
    }
    else
    {
      return false;
    }
  }
  return true;
}
  template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
  bool is_positive_unate_in_i( const TT& tt, uint8_t i )
  {
    auto const tt0 = cofactor0( tt, i );
    auto const tt1 = cofactor1( tt, i );
    for ( auto bit = 0; bit < ( 2 << ( tt.num_vars() - 1 ) ); bit++ )
    {
      if ( get_bit( tt0, bit ) <= get_bit( tt1, bit ) )
      {
        continue;
      }
      else
      {
        return false;
      }
    }
   return true;
  }
  
}







