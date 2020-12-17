#pragma once

#include <stdio.h>
#include <stdlib.h>

/*ILP_SOLVE LIBRARY*/
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
/*Include SOP representation of tt*/
#include "isop.hpp"
#include "traits.hpp"
#include "operations.hpp"
#include "properties.hpp"
#include "implicant.hpp"
#include "print.hpp"
#include "cube.hpp"
#include "constructors.hpp"
#include "bit_operations.hpp"

namespace kitty {

    template<typename TT, typename = std::enable_if_t <kitty::is_complete_truth_table<TT>::value>>
    bool is_negative_unate_in_i(const TT &tt, uint8_t i) {
        auto const tt0 = cofactor0(tt, i);
        auto const tt1 = cofactor1(tt, i);
        for (auto bit = 0; bit < (2 << (tt.num_vars() - 1)); bit++) {
            if (get_bit(tt0, bit) >= get_bit(tt1, bit)) {
                continue;
            } else {
                return false;
            }
        }
        return true;
    }

    template<typename TT, typename = std::enable_if_t <kitty::is_complete_truth_table<TT>::value>>
    bool is_positive_unate_in_i(const TT &tt, uint8_t i) {
        auto const tt0 = cofactor0(tt, i);
        auto const tt1 = cofactor1(tt, i);
        for (auto bit = 0; bit < (2 << (tt.num_vars() - 1)); bit++) {
            if (get_bit(tt0, bit) <= get_bit(tt1, bit)) {
                continue;
            } else {
                return false;
            }
        }
        return true;
    }


    template<typename TT, typename = std::enable_if_t <kitty::is_complete_truth_table<TT>::value>>

    bool is_threshold(const TT &tt, std::vector <int64_t> *plf = nullptr) {

        std::vector <int64_t> linear_form;
        std::vector <int8_t> flipp;  //to store variables where the function is negative unate
        TT tt0 = tt;
        for (auto i = 0u; i < tt.num_vars(); i++) {

            
            if (is_positive_unate_in_i(tt, i) + is_negative_unate_in_i(tt, i) == 0) {
                return false;
            }
            if (is_negative_unate_in_i(tt, i)==1) {
                flip_inplace(tt0, i);
                flipp.push_back(i);
            }

        }


       
        
	
	lprec *lp;
        int Ncol, *colno = NULL, ret = 0;
        REAL *row = NULL;
        Ncol = tt.num_vars() + 1;
	lp = make_lp(0, Ncol);
        if (lp == NULL) {
            return false; 
        }

        colno = (int *) malloc(Ncol * sizeof(*colno));
        row = (REAL *) malloc(Ncol * sizeof(*row));

        if ((colno == NULL) || (row == NULL)) {
            return false;
        }

        if (ret == 0) {

            for (auto j = 0; j <= tt.num_vars(); j++) {
		int k=0;
                colno[k] = j + 1;
                row[j] = 1;
                add_constraintex(lp, k+1, row, colno, GE, 0); 
            }
        }

        if (ret == 0) {
	    std::vector <kitty::cube> pos_const = isop(tt0);
            for (auto k = 0; k < pos_const.size(); k++) {

                kitty::cube pos = pos_const.at(k);

                for (auto i = 0; i < tt.num_vars(); i++) {
                colno[i] = i + 1;
		if (pos.get_mask(i) && pos.get_bit(i))//If the variable exists and it is true
 		   {
                        row[i] = 1;
                    } else {
                        row[i] = 0;
                    }

                }

                colno[tt.num_vars()] = tt.num_vars() + 1; 
                row[tt.num_vars()] = -1;
		add_constraintex(lp, Ncol, row, colno, GE, 0);

            }

        }


        if (ret == 0) {
	    std::vector <kitty::cube> neg_const = isop(unary_not(tt0));
            for (auto k = 0; k < neg_const.size(); k++) {
		kitty::cube neg = neg_const.at(k);
                for (auto i = 0; i < tt.num_vars(); i++) {
                    colno[i] = i + 1;
                     
		if (!neg.get_mask(i) || (neg.get_mask(i) && neg.get_bit(i)))//if the variable is not in the cube or if it is in the cube and it is not in complementary form
		   {
                        row[i] = 1;
                    }
		 else {
                        row[i] = 0;
                    }
                }
                colno[tt.num_vars()] = tt.num_vars() + 1; 
                row[tt.num_vars()] = -1;


                add_constraintex(lp, Ncol, row, colno, LE, -1);

            }

        }

        if (ret == 0) {
            set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
            for (int j = 0; j < Ncol; j++) {
                colno[j] = j + 1;
                row[j] = 1;
            }


            set_obj_fnex(lp, Ncol, row, colno);

        }


        
        set_minim(lp);
	set_verbose(lp, IMPORTANT);
	ret = solve(lp);

        if (ret != 0) {
		return false;
        } else {

	get_variables(lp, row);
	for (auto i = 0; i < Ncol; i++) {
                linear_form.push_back(row[i]); 
            }

            
        for (auto i : flipp) {
                linear_form[i] = -linear_form[i];
                linear_form[Ncol - 1] = linear_form[Ncol - 1] + linear_form[i];
            }


        
        if (row != NULL){
                free(row);
	}
         if (colno != NULL){
                free(colno);
	}
            if (lp != NULL){
		/* clean up such that all used memory by lpsolve is freed */
                delete_lp(lp);
	}

	*plf = linear_form;
        return true;
        }
    }

}

