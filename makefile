# Makefile for meson_deca.

# This is my first makefile and it looks horrible. My bad! =)

##
# Install the package.
##
install:
	# Tell the STAN makefile to use cmdstan (meson_deca/..) 
	# folder for #include command.
	sed -i "s@-Wall@-isystem $\(STANAPI_HOME\).. &@" ../makefile; \
	make reload_libraries


##
# Reload meson_deca libraries into STAN (heavily dependent on CmdStan version.)
##
reload_libraries:
	# This pattern marks the end of the lines we add to the STAN files
	# Note: if you change this, you have to manually delete the lines
	# containing the old pattern.
	EOL_MARK="  // MDECA_LIB"; \
        #
	# Add library functions to the STAN 'functions.hpp' file.
	# Delete all lines containing EOL_MARK
	sed -ie "\@$EOL_MARK@d" ../stan/src/stan/math/prim/mat.hpp; \
	# Insert the info at the end of the file (before '#endif')
	sed -i "s@#endif@#include <meson_deca/lib/c_lib/stan_callable/complex_callable.hpp>$EOL_MARK\n#include <meson_deca/lib/c_lib/model.hpp>$EOL_MARK\n&@" ../stan/src/stan/math/prim/mat.hpp; \
        #
	# Make the necessary changes in 'gm/function_signatures.h'
	sed -ie "\@$EOL_MARK@d" ../stan/src/stan/lang/function_signatures.h; \
        #
	sed -i "s@primitive_types.push_back(DOUBLE_T);@&\nadd(\"A_c\",expr_type(DOUBLE_T,1U),INT_T,VECTOR_T);$EOL_MARK\nadd(\"A_cv\",expr_type(VECTOR_T,1U),VECTOR_T);$EOL_MARK\nadd(\"c_one\",expr_type(DOUBLE_T,1U),DOUBLE_T);$EOL_MARK\nadd(\"c_complex\",expr_type(DOUBLE_T,1U),DOUBLE_T, DOUBLE_T);$EOL_MARK\nadd(\"c_mult\",expr_type(DOUBLE_T,1U),expr_type(DOUBLE_T,1U),expr_type(DOUBLE_T,1U));$EOL_MARK\nadd(\"c_sq_mag\",DOUBLE_T,expr_type(DOUBLE_T,1U));$EOL_MARK\nadd(\"cv_mult\",expr_type(VECTOR_T,1U),expr_type(VECTOR_T,1U),expr_type(VECTOR_T,1U));$EOL_MARK\nadd(\"f_model\",DOUBLE_T,expr_type(VECTOR_T,1U),expr_type(VECTOR_T,1U));$EOL_MARK\nadd(\"cv_sum\",expr_type(DOUBLE_T,1U),expr_type(VECTOR_T,1U));$EOL_MARK\nadd(\"Norm\",DOUBLE_T,expr_type(VECTOR_T,1U),expr_type(MATRIX_T,1U));$EOL_MARK\nadd(\"num_resonances\",INT_T);$EOL_MARK\nadd(\"num_variables\",INT_T);$EOL_MARK@" ../stan/src/stan/lang/function_signatures.h; \
        #
	# STAN binaries must be rebuild
	cd ..;          \
	make clean-all; \
	cd meson_deca
