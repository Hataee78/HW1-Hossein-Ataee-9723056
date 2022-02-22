#include "hw1.h"
#include<iostream>
#include <vector>
#include<random>
#include<functional>
#include <stdexcept>
#include<iomanip>
// namespace algebra{
    Matrix algebra::zeros(std::size_t n, std::size_t m){
        
        Matrix vec ( n , std::vector<double> (m, 0));
        return vec;
    }
    Matrix algebra::ones(std::size_t n, std::size_t m){
        Matrix vec ( n , std::vector<double> (m, 1));
        return vec;
    }
    Matrix algebra::random(size_t n, size_t m, double min, double max){
        if (max <= min)
            throw std::invalid_argument( "logic error" );
        Matrix vec ( n , std::vector<double> (m, 0));
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(min,max);
        auto dice = std::bind ( distribution, generator );
        for(int i=0;i<n;i++)
            for(int j=0;j<m;j++)
                vec[i][j]= dice();
            
        return vec;
    }
    void algebra::show(const Matrix& matrix){
        for(int i=0;i != matrix.size();i++)
            for(int j=0;j != matrix[i].size();j++)
                std::cout<<std::setprecision(3)<< matrix[i][j];
    }
    Matrix algebra:: multiply(const Matrix& matrix, double c){
        Matrix a;
        a = matrix;
        for(int i=0;i != matrix.size();i++)
           for(int j=0;j != matrix[i].size();j++)
              a[i][j]=c*a[i][j] ;
              return a ;    
    }
    Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2){
        int m,n1,n2,h{};
        Matrix e1{};
        m = matrix1.size() ;
        n2 = matrix2.size() ;
        if (m==0||n2==0)
            return e1;
        h =matrix2[0].size(); 
        n1 = matrix1[0].size();
        
        Matrix a;
        Matrix b ( m , std::vector<double> (h, 1));
       
        
        if (n1!=n2)
            throw std::logic_error( "dimensions are wrong" );
        if (m==0&&n1==0&&n2==0&&h==0)    
            return e1;
        if (m==0||n1==0||h==0)  
            throw std::logic_error( "matrix is empty" );

        for(int i=0; i<m;i++)
    {
            for(int j=0;j<h;j++)
        {
                  b[i][j] = 0;

                    for(int k=0;k<n1;k++)
            {
                        b[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
            return b;
    }
    Matrix algebra::sum(const Matrix& matrix, double c){
        
        Matrix a;
        int m = matrix.size() ;
        if (m==0)
            return matrix;
        int n = matrix[0].size();
        if (m==0 || n==0)
            //throw std::logic_error( "matrix is empty" );
            return matrix;
            
            
        a = matrix;
        for(int i=0;i != matrix.size();i++)
           for(int j=0;j != matrix[i].size();j++)
              a[i][j]=c+a[i][j] ;
              return a ;

    }
    Matrix algebra:: sum(const Matrix& matrix1, const Matrix& matrix2){
        int m,n1,n2,h{};
        Matrix e1{};
        
        m = matrix1.size() ;
        n2 = matrix2.size() ;
        if (m!=n2)
            throw std::logic_error( "dimensions are wrong" );
        if (m==0||n2==0)
            return matrix1;
        h =matrix2[0].size(); 
        n1 = matrix1[0].size();
        
        if (h!=n1)
            throw std::logic_error( "dimensions are wrong" );    
        if (m==0&&n1==0&&n2==0&&h==0)    
            return matrix1;
        if (m==0||n1==0||h==0)  
            throw std::logic_error( "matrix is empty" );
            
         Matrix a ( m , std::vector<double> (n1, 1));   
        for(int i=0;i != matrix1.size();i++)
           for(int j=0;j != matrix1[i].size();j++)
              a[i][j]=matrix1[i][j]+matrix2[i][j] ;
              return a ;    
    }
    Matrix algebra:: transpose(const Matrix& matrix){
        int m = matrix.size() ;
        if (m==0)
            return matrix;
        int n = matrix[0].size();
        if (m==0 || n==0)
            return matrix;
            
        
        Matrix a ( n , std::vector<double> (m, 1)) ;       
        
        for(int i=0;i < matrix.size();i++)
           for(int j=0;j < matrix[i].size();j++)
              a[j][i]=matrix[i][j] ;

      return a;        
    }
Matrix algebra::minor(const Matrix& matrix, std::size_t n, std::size_t m){ 
    int minor_row, minor_col{};
    int matrix_col{};
    int matrix_row = matrix.size();
    if (matrix_row > 0)
        matrix_col = matrix[0].size();
    else
        matrix_col = 0;
      
    Matrix a ( (matrix_row-1) , std::vector<double> (matrix_col-1, 1)) ;
        
    for (int i = 0; i < matrix_row; i++) {
        minor_row = i;
        if (i>n)
            minor_row--;
        for (int j = 0; j < matrix_col; j++) {
            minor_col = j;
            if (j>m)
                minor_col--;
            if (i != n && j != m)
                a[minor_row][minor_col] = matrix[i][j];
        }
    }
    return a;
}
    double algebra:: determinant(const Matrix& matrix){
    double det,minor_det{};
    int matrix_col{};
    int matrix_row = matrix.size();
    if (matrix_row > 0)
        matrix_col = matrix[0].size();
    else
        matrix_col = 0;
        
    if (matrix_row!=matrix_col)  
        throw std::logic_error( "matrix is not square" ); 
    if (matrix_col==0)
        return 1;    
        
    if (matrix_col==2)
        return ((matrix[1][1]*matrix[0][0])-(matrix[1][0]*matrix[0][1]));
    else{
        det=0;   
        for(int i=0;i<matrix.size();i++){
            minor_det = matrix[0][i]*algebra::determinant(algebra::minor(matrix,0,i));
            
            if (i%2==0)
                det+=minor_det;
            else 
                det-=minor_det ;   

        }
        
         return det;
       }      
              
    }    
    


    
    Matrix algebra:: inverse(const Matrix& matrix){

    }
     Matrix algebra::concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis){ 
    int matrix_col1{}; 
    int matrix_row1 = matrix1.size(); 
    if (matrix_row1 > 0) 
        matrix_col1 = matrix1[0].size(); 
    else 
        matrix_col1 = 0;  
 
    int matrix_col2{}; 
    int matrix_row2 = matrix2.size(); 
    if (matrix_row2 > 0) 
        matrix_col2 = matrix2[0].size(); 
    else 
        matrix_col2 = 0;   
        if (axis==0) 
            if (matrix_col1!=matrix_col2) 
                throw std::logic_error( "cant concatenate" ); 
        if (axis==1) 
            if (matrix_row1!=matrix_row2) 
                throw std::logic_error( "cant concatenate" );            
     
    if (axis ==1){ 
        Matrix a ( matrix_row1 , std::vector<double> (matrix_col1+matrix_col2, 1)) ; 
        for (int i = 0; i < matrix_row1; i++) { 
        for (int j = 0; j < matrix_col1; j++) { 
  
            // To store elements 
            // of matrix 1 
            a[i][j] = matrix1[i][j]; 
  
            // To store elements 
            // of matrix 2 
            a[i][j + matrix_col1] = matrix2[i][j]; 
            } 
        } 
        return a; 
    }    
         
    if (axis ==0){  
        Matrix a ( matrix_row1+matrix_row2 , std::vector<double> (matrix_col1, 1)) ;  
        for (int i = 0; i < matrix_row1 + matrix_row2; i++) {  
            for (int j = 0; j < matrix_col1; j++) {  
    
                // To store elements  
                // of matrix 1  
                if(i < matrix_row1)
                    a[i][j] = matrix1[i][j];  
                else
                    a[i][j] = matrix2[i-matrix_row1][j];  
            }      
    }
    return a; 
    }            
 
     }
     Matrix algebra::ero_swap(const Matrix& matrix, std:: size_t r1, std:: size_t r2){

        int matrix_col;
        int matrix_row = matrix.size(); 
        if (matrix_row > 0) 
            matrix_col = matrix[0].size(); 
        else 
            return matrix;

        if (matrix_row-1<r1 || matrix_row-1<r2)
            throw std::logic_error( "wrong input" );
            
        Matrix a ( matrix_row , std::vector<double> (matrix_col, 1)) ;
        
           
        for (int i = 0; i < matrix_row; i++) {
            if (i!=r1 && i!=r2){
                a[i]=matrix[i];   
            }
            else{
                a[r1]=matrix[r2];
                a[r2]=matrix[r1];
            }

        }
        return a;  
     }
Matrix algebra:: ero_multiply(const Matrix& matrix, size_t r, double c){
    int matrix_col;
        int matrix_row = matrix.size(); 
        if (matrix_row > 0) 
            matrix_col = matrix[0].size(); 
        else 
            return matrix;

        if (matrix_row-1<r )
            throw std::logic_error( "wrong input" );
            
        Matrix a ( matrix_row , std::vector<double> (matrix_col, 1)) ;
        
           
        for (int i = 0; i < matrix_row; i++) {
            if (i!=r){
                a[i]=matrix[i];   
            }
            else{
                for(int j=0;j<matrix_col;j++){
                    a[i][j]=matrix[i][j]*c;
                }

                
            }

        }
        return a;  
     }