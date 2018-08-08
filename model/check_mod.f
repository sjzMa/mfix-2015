

      MODULE check


        Use param
        Use param1

!                        variables for check_mass_balance
        DOUBLE PRECISION start_time, report_time

        DOUBLE PRECISION accumulation_g, accumulation_X_g(DIM_N_g)
        DOUBLE PRECISION Integral_SUM_R_g, Integral_R_g(DIM_N_g)
        DOUBLE PRECISION flux_out_g(DIMENSION_BC), flux_in_g(DIMENSION_BC)
        DOUBLE PRECISION flux_out_X_g(DIMENSION_BC,DIM_N_g), flux_in_X_g(DIMENSION_BC,DIM_N_g)

        DOUBLE PRECISION accumulation_s(DIM_M), accumulation_X_s(DIM_M, DIM_N_s)
        DOUBLE PRECISION Integral_SUM_R_s(DIM_M), Integral_R_s(DIM_M, DIM_N_s)
        DOUBLE PRECISION flux_out_s(DIMENSION_BC, DIM_M), flux_in_s(DIMENSION_BC, DIM_M)
        DOUBLE PRECISION flux_out_X_s(DIMENSION_BC,DIM_M, DIM_N_s), flux_in_X_s(DIMENSION_BC,DIM_M, DIM_N_s)

      END MODULE check
