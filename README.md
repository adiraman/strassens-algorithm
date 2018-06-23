## Strassen's Subcubic Matrix Multiplication Algorithm

Input :- 2 n x n square matrices

Output :- Compute their matrix multiplication

    | a | b | e | f |
    | c | d | g | h |
	
	  | ae + bg | af + bh | P5 + P4 - P2 - P6 | P1 + P2           |
      | ce + dg | cf + dh | P3 + P4           | P1 + P5 - P3 - P7 |

Algorithm :- 

1. Recursively compute 7 cleverly chosen products.
2. P1 = a(f - h); P2 = (a + b)h; P3 = (c + d)e; P4 = d(g - e); P5 = (a + d)(e + h); P6 = (b - d)(g + h); P7 = (a - c)(e + f).
3. Do necessary additions.
4. For a 2 x 2 or 3 x 3 matrix base case, use brute force method. 
