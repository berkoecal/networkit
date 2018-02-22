#ifndef SYMMATRIX_H
#define SYMMATRIX_H
#include<iostream>
#include<cassert>
#include<vector>
#include<cmath>
#include <stdexcept>
  
  template<typename T, typename IndexType = size_t>
  class SymMatrix{
  public:
    
      //constructors:
    
      /// @description Creates an n times n matrix
      SymMatrix(IndexType n);
      /// @description Creates a n times n matrix with entries initial_value
      SymMatrix(IndexType n, T initial_value );
            
      //operators:
      
      /// @description call operator
      T& operator()(IndexType row,IndexType column);
      
      /// @description increment operator
      void operator+=(const SymMatrix<T, IndexType> & A);
      
      /// @description decrement operator
      void operator-=(const SymMatrix<T, IndexType> & A);
      
      /// @description i dont know the name of this operator
      void operator*=(const T & lambda);
      
      /// @description plus operator
      SymMatrix<T, IndexType> operator+( const SymMatrix<T, IndexType>& A ) const;
      
      /// @description minus operator
      SymMatrix<T, IndexType> operator-( const SymMatrix<T, IndexType>& A ) const;
            
      /// @description multiplication operator
      SymMatrix<T, IndexType> operator*( T lambda ) const;
      
      // Some other elementar stuff:
      
      /// @description get the value of a specified element
      T element(IndexType row, IndexType column) const;
	  
	  ///  @description Additional set function. The access operator() doesn't work for bools
	  void set(IndexType row, IndexType column, T value);
      
      /// @description get the number of columns respective the number of rows
      IndexType size()const{ return _n;}

      //Other Stuff:
      
      /// @return Copy of matrix
      SymMatrix<T, IndexType> copy()const {return SymMatrix<T, IndexType>(*this);}
      
  private:
      IndexType _n;
      std::vector<T> _content;
  };
  
  ///@description one more multiplication
  template<typename T, typename IndexType>
      SymMatrix<T, IndexType> operator*(T lambda, SymMatrix<T, IndexType> A);
  
  /// @description Matrix to ostream
  template<typename T, typename IndexType>
  std::ostream & operator<<(std::ostream& out,const SymMatrix<T, IndexType> & A);
 
  
 
  /*****************************************************************************************************************************************************************************************************
   * 											Code Section
   *****************************************************************************************************************************************************************************************************/


	  template<typename T, typename IndexType>
      SymMatrix<T, IndexType>::SymMatrix(IndexType n):_n(n),_content(static_cast<IndexType>(n*(n+1)/2),T(0.0))
      {
		if(n==0)
			throw std::range_error("Not allowed! Try to construct an empty matrix.");
      }
      template<typename T, typename IndexType>
      SymMatrix<T, IndexType>::SymMatrix(IndexType n,T initial_value ):_n(n),_content(static_cast<IndexType>((n*(n+1))/2),initial_value)
      {
		if(n==0)
			throw std::range_error("Unallowed try to construct an empty matrix.");
      }
   
      template<typename T, typename IndexType>
      T& SymMatrix<T, IndexType>::operator()(IndexType row,IndexType column)
      {
		assert(row<_n and column<_n);
		if (row > column)
			std::swap(row,column);
		return _content[static_cast<IndexType>(column*(column+1)/2)+row];
      }
      
      template<typename T, typename IndexType>
      T SymMatrix<T, IndexType>::element(IndexType row, IndexType column) const
      {
		assert(row<_n and column<_n);
		if (row > column)
			std::swap(row,column);
		return _content[static_cast<IndexType>((column*(column+1))/2)+row];
      }
	  
	  template<typename T, typename IndexType>
	  void SymMatrix<T, IndexType>::set(IndexType row, IndexType column, T new_value)
	  {
		assert(row<_n and column<_n);
		if (row > column)
			std::swap(row,column);
		_content[static_cast<IndexType>((column*(column+1))/2)+row] = new_value;
	  }
  
      template<typename T, typename IndexType>
      void SymMatrix<T, IndexType>::operator+=(const SymMatrix<T, IndexType> & A)
	  {
		assert(_n == A._n);
		for(IndexType k = 0; k < _content.size();++k)
			_content[k]+=A._content[k];
	  }
      
      template<typename T, typename IndexType>
      void SymMatrix<T, IndexType>::operator-=(const SymMatrix<T, IndexType> & A)
	  {
		assert(_n == A._n);
		for(IndexType k = 0; k < _content.size();++k)
			_content[k]-=A._content[k];
      }
      
      template<typename T, typename IndexType>
      void SymMatrix<T, IndexType>::operator*=(const T & lambda){
	  for(IndexType k = 0;k<_content.size();++k)
	      _content[k]*=lambda;
      }
      
      template<typename T, typename IndexType> 
      SymMatrix<T, IndexType> SymMatrix<T, IndexType>::operator+( const SymMatrix<T, IndexType>& A ) const
      {
		assert(_n==A._n);
		SymMatrix<T, IndexType> Res(_n);
		for(IndexType k = 0;k<_content.size();++k)
			Res._content[k]=_content[k]+A._content[k];
		return Res;
      }
      
      template<typename T, typename IndexType> 
      SymMatrix<T, IndexType> SymMatrix<T, IndexType>::operator-( const SymMatrix<T, IndexType>& A ) const
      {
		assert(_n==A._n);
		SymMatrix<T, IndexType> Res(_n);
		for(IndexType k = 0;k<_content.size();++k)
			Res._content[k]=_content[k]-A._content[k];
		return Res;
      }
      
      
      template<typename T, typename IndexType>
      SymMatrix<T, IndexType> SymMatrix<T, IndexType>::operator*(T lambda) const
      {
	   SymMatrix<T, IndexType> Res(_n,0);
	   for(IndexType k = 0;k< _content.size();++k)
	   {
	      Res._content[k] = lambda* _content[k];
	   }
	   return Res;
      }
      
      template<typename T, typename IndexType>
      SymMatrix<T, IndexType> operator*(T lambda, SymMatrix<T, IndexType> A)
      {
		return (A*lambda);
      }


#endif

