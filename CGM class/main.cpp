#include <iostream>
#include "distribution_function.h"

using namespace std;
int main() {
	distribution_function_space::distribution_function_D2Q9 f(3, 3);
	f(2, 3, 1, 2);


	return 0;
}