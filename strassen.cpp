#include "strassen.h"

matrix_t add(const matrix_t& a, const matrix_t& b)
{
	int rows = a.size();

	vector_t r(rows);
	matrix_t c(rows,r);

	for(int i=0; i<rows; ++i)
	{
		for(int j=0; j<rows; ++j)
		{
			c[i][j] = a[i][j] + b[i][j];
		}
	}
	return c;
}

matrix_t subtract(const matrix_t& a, const matrix_t& b)
{
	int rows = a.size();

	vector_t r(rows);
	matrix_t c(rows,r);

	for(int i=0; i<rows; ++i)
	{
		for(int j=0; j<rows; ++j)
		{
			c[i][j] = a[i][j] - b[i][j];
		}
	}
	return c;
}

matrix_t multiply(const matrix_t& a, const matrix_t& b)
{
	int rows = a.size();
	vector_t temp(rows);
	matrix_t c(rows, temp);

	for(int i=0; i<rows; ++i)
	{
		for(int k=0; k<rows; ++k)
		{
			for(int j=0; j<rows; ++j)
			{
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return c;
}

bool equals(const matrix_t& a, const matrix_t& b)
{
	int rows = a.size();
	int counter = 0;

	for(int i=0; i<rows; ++i)
	{
		for(int j=0; j<rows; ++j)
		{
			if(a[i][j] == b[i][j]) ++counter;
		}
	}

	return (counter == rows*rows);
}

void print(const matrix_t& a)
{
	int rows;
	rows = a.size();

	for(int i=0; i<rows; ++i)
	{
		for(int j=0; j<rows; ++j)
		{
			std::cout << a[i][j] << ", " << std::flush;
		}
		std::cout << std::endl;
	}
}

matrix_t strassen(const matrix_t& a, const matrix_t& b)
{
	int rows = a.size();

	// Check for base case
	if(rows <= 2)
	{
		return multiply(a, b);
	}

	int newRows = rows/2;
	vector_t inner(newRows);

	// initialize the 8 divided sub matrices, products and sums
	matrix_t a11(newRows, inner), a12(newRows, inner),
			 a21(newRows, inner), a22(newRows, inner),
			 b11(newRows, inner), b12(newRows, inner),
			 b21(newRows, inner), b22(newRows, inner),
			 C11(newRows, inner), C12(newRows, inner),
			 C21(newRows, inner), C22(newRows, inner),
			 aRes(newRows, inner), bRes(newRows, inner),
			 P1(newRows, inner),P2(newRows, inner),P3(newRows, inner),
			 P4(newRows, inner),P5(newRows, inner),P6(newRows, inner),
			 P7(newRows, inner);

	// Divide into sub matrices
	for(int i=0; i<newRows; ++i)
	{
		for(int j=0; j<newRows; ++j)
		{
			a11[i][j] = a[i][j];
			a12[i][j] = a[i][j+newRows];
			a21[i][j] = a[i+newRows][j];
			a22[i][j] = a[i+newRows][j+newRows];

			b11[i][j] = b[i][j];
			b12[i][j] = b[i][j+newRows];
			b21[i][j] = b[i+newRows][j];
			b22[i][j] = b[i+newRows][j+newRows];
		}
	}

	// Calculate 7 Product matrices recursively
	P1 = strassen(add(a11,a22), add(b11,b22));//p1=(a11+a22)*(b11+b22)
	P2 = strassen(add(a21, a22), b11);//p2=(a21+a22)*b11
	P3 = strassen(a11, subtract(b12, b22));//p3=a11*(b12-b22)
	P4 = strassen(a22, subtract(b21, b11));//p4=a22*(b21-b11)
	P5 = strassen(add(a11, a12), b22);//p5=(a11+a12)*b22
	P6 = strassen(subtract(a21,a11), add(b11,b12));//p6=(a21-a11)*(b11+b12)
	P7 = strassen(subtract(a12,a22), add(b21,b22));//p7=(a12-a22)*(b21+b22)

	// Calculate the result sub matrices from p1-p7
	C11 = subtract(add(add(P1, P4), P7), P5);//c11 = p1+p4+p7-p5
	C12 = add(P3, P5);//c12 = p3+p5
	C21 = add(P2, P4);//c21 = p2+p4
	C22 = subtract(add(add(P1, P3), P6), P2);//c22 = p1+p3+p6-p2

	// Assimilate the rhs matrix
	vector_t innerN(rows);
	matrix_t c(rows,innerN);

	for(int i=0; i<newRows; ++i)
	{
		for(int j=0; j<newRows; ++j)
		{
			c[i][j] = C11[i][j];
			c[i][j+newRows] = C12[i][j];
			c[i+newRows][j] = C21[i][j];
			c[i+newRows][j+newRows] = C22[i][j];
		}
	}

	return c;
}

