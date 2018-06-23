/*
 * =====================================================================================
 *
 *       Filename:  strassen.h
 *
 *    Description:  Header file for strassen matrix multiplication
 *
 *        Version:  1.0
 *        Created:  10/21/2017 10:19:35 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Aditya Raman
 *
 * =====================================================================================
 */

#include <cstdlib>
#include <iostream>
#include <vector>

typedef std::vector< std::vector<int> > matrix_t;
typedef std::vector<int> vector_t;

matrix_t add(const matrix_t& a, const matrix_t& b);

matrix_t subtract(const matrix_t& a, const matrix_t& b);

matrix_t multiply(const matrix_t& a, const matrix_t& b);

bool equals(const matrix_t& a, const matrix_t& b);

void print(const matrix_t& a);

matrix_t strassen(const matrix_t& a, const matrix_t& b);

