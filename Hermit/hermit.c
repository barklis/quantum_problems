#include <stdio.h>

// Recursive function for calculation of Hermite polynomial
// Based on equation H(n)[x] = 2 * x * H(n-1)[x] - 2 * (n-1) * H(n-2)[x]
// Should be detection of overflow somehow done...
long long hermite(long long x, long long n){ //Calculate n-th Hermite polynomial of x
	if(n == 0)		
		return 1;	// First two known
	else if(n == 1)		// values of Hermite polynomial
		return 2 * x;	// H(0)[x] = 1 and H(1)[x] = 2x
	else
		return 2 * x * hermite(x, n - 1) - 2 * (n - 1) * hermite(x, n - 2); // Calculate H(n)[x]
}	

//Maybe sometime will be better then previous
int solve_hemit(int *tab, int n){
	return 0;
}

//Tests
int main(){
	long long degree, value;
	printf("Enter degree of Hermit polynomial: ");
	scanf("%lld", &degree);
	printf("Enter argument value of polynomial: ");
	scanf("%lld", &value);
	printf("Value of polynomial is %lld\n", hermite(value, degree));
	return 0;
}
