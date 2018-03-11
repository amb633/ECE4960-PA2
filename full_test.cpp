#include <iostream>
#include <cstdlib>
#include <cmath>

#include "full_func.hpp"
#include "test_full_func.hpp"

using namespace std;
using namespace full;
using namespace test_full;

int main(int argc, char const *argv[])
{
	cout << endl << boolalpha;
	cout << "test for matrix retrieve element : " << test_matrix_retrieve_element() << endl;
	cout << "test for vector retrieve element : " << test_vector_retrieve_element() << endl;
	cout << "test for change element function : " << test_change_element() << endl;
	cout << "test for copy matrix function    : " << test_copy_matrix() << endl;
	cout << "test for scalar multiple function: " << test_scalar_multiple() << endl;
	cout << "test for row permute function    : " << test_row_permute() << endl;
	cout << "test for row scale function      : " << test_row_scale() << endl;
	cout << "test for matrix product function : " << test_matrix_product() << endl;
	cout << "test for norm caluclation        : " << test_calculate_norm() << endl;
	cout << "test for matrix decomposition    : " << test_matrix_decomposition() << endl;
	cout << endl;
	return 0;
}