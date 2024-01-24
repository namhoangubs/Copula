library("R6")
library("assertive")

Copula_Class <- R6Class("Copula_Class",
                        
                        private = list(
                                ..Confidence_Level = 0.05,
                                ..Rotations = FALSE,
                                ..Method = "mle",
                                ..Copula_Family = c(1,2,3,4,5,6),
                                ..Selection_Criteria = "AIC"
                        ),
                        
                        active = list(
                                
                                Confidence_Level = function(value){
                                        if(missing(value)){
                                                private$..Confidence_Level
                                        } else{
                                                assert_is_numeric(value)
                                                assert_all_are_in_closed_range(value, lower = 0, upper = 1)
                                                private$..Confidence_Level <- value
                                        }
                                        
                                },
                                
                                Rotations = function(value){
                                        if(missing(value)){
                                                private$..Rotations
                                        } else{
                                                assert_is_a_bool(value)
                                                private$..Rotations <- value
                                        }
                                },
                                
                                Method = function(value){
                                        if(missing(value)){
                                                private$..Method
                                        } else{
                                                assert_is_a_string(value)
                                                private$..Method <- value
                                        } 
                                        
                                },
                                
                                Copula_Family = function(value){
                                        if(missing(value)){
                                                private$..Copula_Family
                                        } else{
                                                assert(assert_is_vector(value), assert_is_numeric(value), combine = "or")
                                                private$..Copula_Family <- value
                                        } 
                                },
                                
                                Selection_Criteria = function(value){
                                        if(missing(value)){
                                                private$..Selection_Criteria
                                        } else{
                                                assert_is_a_string(value)
                                                private$..Selection_Criteria <- value
                                        }
                                }
                                
                        ),
                        
                        public = list(
                                
                                Copula_Calibration = function(index_input, first_series, second_series){
                                        u1 <- pobs(index_input[,first_series])
                                        u2 <- pobs(index_input[,second_series])
                                        return (BiCopSelect(u1,u2,selectioncrit = private$..Selection_Criteria,
                                                                    familyset = private$..Copula_Family, level = private$..Confidence_Level,
                                                                    rotations = private$..Rotations, method =private$..Method))
                                },
                                
                                Copula_Generator = function(Cop_family, Cop_par, Cop_par2, no_of_simulation){
                                        Simulated_U_Value = BiCopSim(no_of_simulation, family = Cop_family,
                                                                     par = Cop_par, par2 = Cop_par2)
                                        return (BiCopHinv(Simulated_U_Value[,1], Simulated_U_Value[,2], 
                                                                  family = Cop_family, par = Cop_par,
                                                          par2 = Cop_par2)
                                                )
                                },
                                
                                Set_Parameters = function(confidence_interval, rotations, method, copula_family){
                                        private$..Confidence_Level <- confidence_interval
                                        private$..Rotations <- rotations
                                        private$..Method <- method
                                        private$..Copula_Family <- copula_family
                                },
                                
                                
                                initialize = function(Confidence_Level, Rotations, Method, Copula_Family){
                                        
                                        if(!missing(Confidence_Level)){
                                                private$..Confidence_Level <- Confidence_Level
                                        }
                                        
                                        if(!missing(Rotations)){
                                                private$..Rotations <- Rotations
                                        }
                                        
                                        if(!missing(Method)){
                                                private$..Method <- Method
                                        }
                                        
                                        if(!missing(Copula_Family)){
                                                private$..Copula_Family <- Copula_Family
                                        }
                                        
                                }
                                
                                
                        )
                        
)

Marginal_Dist_Class <- R6Class("Marginal_Dist_Class",
                          
                               private = list(
                                  ..distribution_name = "norm",
                                  ..df = 10000
                                  
                          ),
                          
                                active = list(
                                        distribution_name = function(value){
                                                if(missing(value)){
                                                  private$..distribution_name
                                                } else{
                                                  assert_is_a_string(value)
                                                  private$..Confidence_Level <- value
                                                }
                                        },
                                  
                                        df = function(value){
                                                if(missing(value)){
                                                        private$..distribution_name
                                                } else{
                                                        assert_is_integer(value)
                                                        private$..df <- value
                                                }
                                        }
                          
                                ),
                          
                                public = list(
                                        norm_dist_generation = function(series, CDF_value){
                                                norm_dist_fitting <- fitdist(series, distr = "norm")
                                                norm_dist_generator <- qnorm(p = CDF_value, mean = norm_dist_fitting$estimate['mean'],
                                                                            sd = norm_dist_fitting$estimate['sd'])
                                                
                                        },
                                
                                        initialize = function(distribution_name, df){
                                                
                                                if(!missing(distribution_name)){
                                                        private$..distribution_name <- distribution_name
                                                }
                                                
                                                if(!missing(df)){
                                                        private$..df<- df
                                                }
                                        }
                                        
                                        )
                          
                                
)

                          
