#include "strassen.h"

int main()
{
	int sz = 128;
	vector_t r1(sz);
	matrix_t a(sz, r1);
	matrix_t b(sz, r1);

	for(int i=0; i<sz; ++i)
	{
		for(int j=0; j<sz; ++j)
		{
			a[i][j] = (rand()%10 + 1);
			b[i][j] = (rand()%10 + 1);
		}
	}

	matrix_t c = strassen(a, b);
	matrix_t c1 = multiply(a, b);
	bool res = equals(c, c1);

	if(res)
	{
		std::cout << "Strassen matches brute force method !" << std::endl;
	}
	else
	{
		std::cout << "Implementation error in Strassen" << std::endl;
	}

	return 0;
}


