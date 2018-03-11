#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "compressed_func.hpp"

using namespace std;

namespace large_matrix_test{
	bool test_row_permute( compressed::comp_r_mat* AC );
	bool test_row_scale( compressed::comp_r_mat* AC );
	bool test_matrix_product( compressed::comp_r_mat* AC );
	bool test_calculate_norm( compressed::comp_r_mat* AC );
	bool test_matrix_decomposition(compressed::comp_r_mat* AC );

	int run_large_matrix_tests();
}