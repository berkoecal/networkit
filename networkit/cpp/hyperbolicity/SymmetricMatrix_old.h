#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>
#include<cassert>
#include<cmath>
#include <stdexcept>
#include <vector>
  
  template<typename T>
  class SymMatrix{
  public:
    
      //constructors:
    
      /// @description Creates an n times n matrix
      SymMatrix(size_t n);
      /// @description Creates a n times n matrix with entries initial_value
      SymMatrix(size_t n, T initial_value );
            
      //operators:
      
      /// @description call operator
      T& operator()(size_t row,size_t column);
      
      /// @description increment operator
      void operator+=(const SymMatrix<T> & A);
      
      /// @description decrement operator
      void operator-=(const SymMatrix<T> & A);
      
      /// @description i dont know the name of this operator
      void operator*=(const T & lambda);
      
      /// @description plus operator
      SymMatrix<T> operator+( const SymMatrix<T>& A ) const;
      
      /// @description minus operator
      SymMatrix<T> operator-( const SymMatrix<T>& A ) const;
            
      /// @description multiplication operator
      SymMatrix<T> operator*( T lambda ) const;
      
      // Some other elementar stuff:
      
      /// @description get the value of a specified element
      T element(size_t row, size_t column) const;
      
      /// @description get the number of columns respective the number of rows
      size_t size()const{ return _n;}
      
      //Other Stuff:
      
      /// @return Copy of matrix
      SymMatrix<T> copy()const {return SymMatrix<T>(*this);}
      
  private:
      size_t _n;
      std::vector<T> _content;
  };
  
  ///@description one more multiplication
  template<typename T>
      SymMatrix<T> operator*(T lambda, SymMatrix<T> A);
  
  /// @description Matrix to ostream
  template<typename T>
  std::ostream & operator<<(std::ostream& out,const SymMatrix<T> & A);
 
  
 
  /*****************************************************************************************************************************************************************************************************
   * 											Code Section
   *****************************************************************************************************************************************************************************************************/


	  template<typename T>
      SymMatrix<T>::SymMatrix(size_t n):_n(n),_content(static_cast<size_t>(n*(n+1)/2),T(0.0))
      {
		if(n==0)
			throw std::range_error("Not allowed! Try to construct an empty matrix.");
      }
      template<typename T>
      SymMatrix<T>::SymMatrix(size_t n,T initial_value ):_n(n),_content(static_cast<size_t>(n*(n+1)/2),initial_value)
      {
		if(n==0)
			throw std::range_error("Unallowed try to construct an empty matrix.");
      }
   
      template<typename T>
      T& SymMatrix<T>::operator()(size_t row,size_t column)
      {
		assert(row<_n and column<_n);
		if (row > column)
			std::swap(row,column);
		return _content[static_cast<size_t>(row*(row+1)/2)+column];
      }
      
      template<typename T>
      T SymMatrix<T>::element(size_t row, size_t column) const
      {
		assert(row<_n and column<_n);
		if (row > column)
			std::swap(row,column);
		return _content[static_cast<size_t>(row*(row+1)/2)+column];
      }
  
      template<typename T>
      void SymMatrix<T>::operator+=(const SymMatrix<T> & A)
	  {
		assert(_n == A._n);
		for(size_t k = 0; k < _content.size();++k)
			_content[k]+=A._content[k];
	  }
      
      template<typename T>
      void SymMatrix<T>::operator-=(const SymMatrix<T> & A)
	  {
		assert(_n == A._n);
		for(size_t k = 0; k < _content.size();++k)
			_content[k]-=A._content[k];
      }
      
      template<typename T>
      void SymMatrix<T>::operator*=(const T & lambda){
	  for(size_t k = 0;k<_content.size();++k)
	      _content[k]*=lambda;
      }
      
      template<typename T> 
      SymMatrix<T> SymMatrix<T>::operator+( const SymMatrix<T>& A ) const
      {
		assert(_n==A._n);
		SymMatrix<T> Res(_n);
		for(size_t k = 0;k<_content.size();++k)
			Res._content[k]=_content[k]+A._content[k];
		return Res;
      }
      
      template<typename T> 
      SymMatrix<T> SymMatrix<T>::operator-( const SymMatrix<T>& A ) const
      {
		assert(_n==A._n);
		SymMatrix<T> Res(_n);
		for(size_t k = 0;k<_content.size();++k)
			Res._content[k]=_content[k]-A._content[k];
		return Res;
      }
      
      
      template<typename T>
      SymMatrix<T> SymMatrix<T>::operator*(T lambda) const
      {
	   SymMatrix<T> Res(_n,0);
	   for(size_t k = 0;k< _content.size();++k)
	   {
	      Res._content[k] = lambda* _content[k];
	   }
	   return Res;
      }
      
      template<typename T>
      SymMatrix<T> operator*(T lambda, SymMatrix<T> A)
      {
		return (A*lambda);
      }


#endif

