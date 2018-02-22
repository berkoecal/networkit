#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>
#include<cassert>
#include<cmath>
#include <stdexcept>

  
  template<typename T>
  class Matrix{
  public:
    
      //constructors:
    
      /// @description Creates an n times n matrix
      Matrix(size_t n);
      /// @description Creates a num_rows times num_columns matrix
      Matrix(size_t num_rows,size_t num_columns );
      /// @description Creates a num_rows times num_columns matrix with entries initial_value
      Matrix(size_t num_rows,size_t num_columns,T initial_value );
            
      //operators:
      
      /// @description call operator
      T& operator()(size_t row,size_t column);
      
      /// @description increment operator
      void operator+=(const Matrix<T> & A);
      
      /// @description decrement operator
      void operator-=(const Matrix<T> & A);
      
      /// @description I don't know the name of this operator
      void operator*=(const T & lambda);
      
      /// @description plus operator
      Matrix<T> operator+( const Matrix<T>& A ) const;
      
      /// @description minus operator
      Matrix<T> operator-( const Matrix<T>& A ) const;
            
      /// @description multiplication operator
      Matrix<T> operator*( T lambda ) const;
      
      // Some other elementary stuff:
      
      /// @description get the value of a specified element
      T element(size_t row, size_t column) const;
      
      // @description set (@param row, @param column)-th entry to @param value
      void set(size_t row, size_t column, T value);

      /// @description get the number of columns
      size_t num_columns()const{ return _num_columns;}
      
      /// @description get the number of columns
      size_t num_rows()const{ return _num_rows;}
      
      /// @return true if the matrix is square
      bool is_square() const{ return (_num_rows == _num_columns);      }
      
      // Elemental matrix operations:
      
      /// @description Switches the rows a and b.
      void switch_rows(size_t a, size_t b);
      
      /// @description multiplies the row't row with the factor lambda
      void multiply_row(size_t row, T lambda);
      
      /// @description add the from-row lambda times to the to row.
      void add_row_to_row(T lambda, size_t from, size_t to);
      
      /// @description Switches the columns a and b.
      void switch_columns(size_t a, size_t b);
      
      /// @description multiplies the column't column with the factor lambda
      void multiply_column(size_t column, T lambda);
      
      /// @description add the from-column lambda times to the to column.
      void add_column_to_column(T lambda, size_t from, size_t to);
      
      //Other Stuff:
      
      /// @return Copy of matrix
      Matrix<T> copy()const {return Matrix<T>(*this);}
      
  private:
      size_t _num_rows;
      size_t _num_columns;
      std::vector<T> _content;
  };
  
  ///@description one more multiplication
  template<typename T>
      Matrix<T> operator*(T lambda, Matrix<T> A);
  
  /// @description Matrix to ostream
  template<typename T>
  std::ostream & operator<<(std::ostream& out,const Matrix<T> & A);
 
  
 
  /*****************************************************************************************************************************************************************************************************
   * 											Code Section
   *****************************************************************************************************************************************************************************************************/
      template<typename T>
      Matrix<T>::Matrix(size_t n):_num_rows(n),_num_columns(n),_content(n*n,T(0.0))
      {
	if(n==0)
	  throw std::range_error("Unallowed try to construct an empty matrix.");
      }

      template<typename T>
      Matrix<T>:: Matrix(size_t num_rows,size_t num_columns ):_num_rows(num_rows),_num_columns(num_columns),_content(num_rows*num_columns,T(0.0))
      {
	  if(num_rows==0 || num_columns==0)
	      throw std::range_error("Unallowed try to construct an empty matrix.");	
      }
      template<typename T>
      Matrix<T>::Matrix(size_t num_rows,size_t num_columns,T initial_value ):_num_rows(num_rows),_num_columns(num_columns),_content(num_rows*num_columns,initial_value)
      {
	  if(num_rows==0 || num_columns==0)
	    throw std::range_error("Unallowed try to construct an empty matrix.");
      }
   
      template<typename T>
      T& Matrix<T>::operator()(size_t row,size_t column)
      {
	assert(row<_num_rows);
	assert(column<_num_columns);
	return _content[row*_num_columns+column];
      }
      
      template<typename T>
      T Matrix<T>::element(size_t row, size_t column) const
      {
	assert(row<_num_rows);
	assert(column<_num_columns);
	return _content[row*_num_columns+column];
      }

	  template<typename T, typename IndexType>
	  void Matrix<T>::set(size_t row, size_t column, T new_value)
	  {
		assert(row<_num_rows and column<_num_columns);
		if (row > column)
			std::swap(row,column);
		_content[static_cast<size_t>(row*_num_columns +column)] = new_value;
	  }
  
      template<typename T>
      void Matrix<T>::operator+=(const Matrix<T> & A){
	  assert(_num_columns == A._num_columns);
	  assert(_num_rows == A._num_rows);
	  for(size_t k = 0;k<_num_rows*_num_columns;++k)
	      _content[k]+=A._content[k];
      }
      
      template<typename T>
      void Matrix<T>::operator-=(const Matrix<T> & A){
	  assert(_num_columns == A._num_columns);
	  assert(_num_rows == A._num_rows);
	  for(size_t k = 0;k<_num_rows*_num_columns;++k)
	      _content[k]-=A._content[k];
      }
      
      template<typename T>
      void Matrix<T>::operator*=(const T & lambda){
	  for(size_t k = 0;k<_num_rows*_num_columns;++k)
	      _content[k]*=lambda;
      }
      
      template<typename T>
      void Matrix<T>::switch_rows(size_t a, size_t b)
      {
	assert(a<_num_rows);
	assert(b<_num_rows);
	for(size_t k = 0;k<_num_columns;++k)
	{
	   T temp =_content[a*_num_columns+k];
	  _content[a*_num_columns+k]=_content[b*_num_columns+k];
	  _content[b*_num_columns+k]=temp;
	}
      }
      
      template<typename T>
      void Matrix<T>::multiply_row(size_t row, T factor)
      {
	  assert(row<_num_rows);
	  
	for(size_t k = 0;k<_num_columns;++k)
	{
	   _content[row*_num_columns+k]*=factor;  
	} 
      }
      
      template<typename T>      
      void Matrix<T>::add_row_to_row(T lambda, size_t from, size_t to)
      {
	assert(from<_num_rows);
	assert(to<_num_rows);
	for(size_t k = 0;k<_num_columns;++k)
	{
	  _content[to*_num_columns+k]+= (lambda* _content[from*_num_columns+k]);
	}
      }
      
      template<typename T>
      void Matrix<T>::switch_columns(size_t a, size_t b)
      {
	assert(a<_num_rows);
	assert(b<_num_rows);
	for(size_t k = 0;k<_num_rows;++k)
	{
	   T temp =_content[k*_num_columns+a];
	  _content[k*_num_columns+a]=_content[k*_num_columns+b];
	  _content[k*_num_columns+b]=temp;
	}
      }
      
      template<typename T>
      void Matrix<T>::multiply_column(size_t column, T factor)
      {
	assert(column<_num_columns);
	  
	for(size_t k = 0;k<_num_rows;++k)
	{
	   _content[column*_num_rows+k]*=factor;  
	} 
      }
      
      template<typename T>      
      void Matrix<T>::add_column_to_column(T lambda, size_t from, size_t to)
      {
	assert(from<_num_columns);
	assert(to<_num_columns);
	for(size_t k = 0;k<_num_rows;++k)
	{
	  _content[k*_num_columns+to]+= (lambda* _content[k*_num_columns+from]);
	}
      }
      
      template<typename T> 
      Matrix<T> Matrix<T>::operator+( const Matrix<T>& A ) const
      {
	 assert(_num_columns==A._num_columns);
	 assert(_num_rows==A._num_rows);
	 Matrix<T> Res(_num_rows,_num_columns);
	 for(size_t k = 0;k<_num_rows*_num_columns;++k)
	    Res._content[k]=_content[k]+A._content[k];
	 return Res;
      }
      
      template<typename T> 
      Matrix<T> Matrix<T>::operator-( const Matrix<T>& A ) const
      {
	 assert(_num_columns==A._num_columns);
	 assert(_num_rows==A._num_rows);
	 Matrix<T> Res(_num_rows,_num_columns);
	 for(size_t k = 0;k<_num_rows*_num_columns;++k)
	    Res._content[k]=_content[k]-A._content[k];
	 return Res;
      }
      
      
      template<typename T>
      Matrix<T> Matrix<T>::operator*(T lambda) const
      {
	   Matrix<T> Res(_num_rows,_num_columns,0);
	   for(size_t k = 0;k<_num_rows*_num_columns;++k)
	   {
	      Res._content[k] = lambda* _content[k];
	   }
	   return Res;
      }
      
      template<typename T>
      Matrix<T> operator*(T lambda, Matrix<T> A)
      {
		return (A*lambda);
      }


#endif

