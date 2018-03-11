#include "full_func.hpp"
#include "compressed_func.hpp"

namespace wilkinson_test{
	
	bool test_retrieve_element( compressed::comp_r_mat* AC , vector< vector<double>>* AF );
	bool test_change_element( compressed::comp_r_mat* AC , vector< vector<double>>* AF );
	bool test_copy_matrix( compressed::comp_r_mat* AC , vector< vector<double>>* AF );
	bool test_scalar_multiple(compressed::comp_r_mat* AC , vector< vector<double>>* AF );

	bool test_row_permute( compressed::comp_r_mat* AC , vector< vector<double>>* AF );
	bool test_row_scale( compressed::comp_r_mat* AC , vector< vector<double>>* AF );
	bool test_matrix_product( compressed::comp_r_mat* AC , vector< vector<double>>* AF , vector<double>* B);
	bool test_calculate_norm( compressed::comp_r_mat* AC , vector< vector<double>>* AF , vector<double>* B);
	bool test_matrix_decomposition( compressed::comp_r_mat* AC , vector< vector<double>>* AF );

	bool test_jacobi_solver( compressed::comp_r_mat* AC , vector< vector<double>>* AF , vector<double>* B);

	int run_wilkinson_tests();
}