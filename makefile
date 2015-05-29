# Makefile for meson_deca.

# This is my first makefile and it looks horrible. My bad! =)

##
# Install the package.
##
install:
	# Tell the STAN makefile to link cmdstan folder 
	# when building C++ files (this allows us to use
        # '#include <meson_deca/..>' statements in C++ files) 
	sed -i "s@-Wall@-isystem $$\(STANAPI_HOME\).. &@" ../makefile; \
	make reload_libraries


##
# Reload meson_deca libraries into STAN (heavily dependent on CmdStan version.
# Current support: CmdStan 2.6.2.)
##
reload_libraries:
	# Add library functions to the STAN 'functions.hpp' file.
	# Delete all lines containing EOL_MARK
	sed -ie "\@  // MDECA_LIB@d" ../stan/src/stan/math/prim/mat.hpp; \
	# Insert the info at the end of the file (before '#endif')
	sed -i "s@#endif@#include <meson_deca/lib/c_lib/stan_callable/complex_callable.hpp>  // MDECA_LIB\n#include <meson_deca/lib/c_lib/model.hpp>  // MDECA_LIB\n&@" ../stan/src/stan/math/prim/mat.hpp; \
        #
	# Make the necessary changes in 'gm/function_signatures.h'
	sed -ie "\@  // MDECA_LIB@d" ../stan/src/stan/lang/function_signatures.h; \
        #
	sed -i "s@primitive_types.push_back(DOUBLE_T);@&\nadd(\"A_c\",expr_type(DOUBLE_T,1U),INT_T,VECTOR_T);  // MDECA_LIB\nadd(\"A_cv\",expr_type(VECTOR_T,1U),VECTOR_T);  // MDECA_LIB\nadd(\"A_v_background_abs2\",VECTOR_T,VECTOR_T);  // MDECA_LIB\nadd(\"c_one\",expr_type(DOUBLE_T,1U),DOUBLE_T);  // MDECA_LIB\nadd(\"c_complex\",expr_type(DOUBLE_T,1U),DOUBLE_T, DOUBLE_T);  // MDECA_LIB\nadd(\"c_mult\",expr_type(DOUBLE_T,1U),expr_type(DOUBLE_T,1U),expr_type(DOUBLE_T,1U));  // MDECA_LIB\nadd(\"c_sq_mag\",DOUBLE_T,expr_type(DOUBLE_T,1U));  // MDECA_LIB\nadd(\"cv_mult\",expr_type(VECTOR_T,1U),expr_type(VECTOR_T,1U),expr_type(VECTOR_T,1U));  // MDECA_LIB\nadd(\"f_model\",DOUBLE_T,expr_type(VECTOR_T,1U),expr_type(VECTOR_T,1U));  // MDECA_LIB\nadd(\"f_model\",DOUBLE_T,expr_type(VECTOR_T,1U),expr_type(VECTOR_T,1U), VECTOR_T, VECTOR_T);  // MDECA_LIB\nadd(\"cv_sum\",expr_type(DOUBLE_T,1U),expr_type(VECTOR_T,1U));  // MDECA_LIB\nadd(\"Norm\",DOUBLE_T,expr_type(VECTOR_T,1U),expr_type(MATRIX_T,1U));  // MDECA_LIB\nadd(\"Norm\",DOUBLE_T,expr_type(VECTOR_T,1U),expr_type(MATRIX_T,1U), VECTOR_T, VECTOR_T);  // MDECA_LIB\nadd(\"num_background\",INT_T);  // MDECA_LIB\nadd(\"num_resonances\",INT_T);  // MDECA_LIB\nadd(\"num_variables\",INT_T);  // MDECA_LIB@" ../stan/src/stan/lang/function_signatures.h; \
        #
	# STAN binaries must be rebuild
	cd ..;          \
	make clean-all; \
	cd meson_deca
