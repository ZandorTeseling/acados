/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#ifndef ACADOS_SIM_SIM_DISCRETE_H_
#define ACADOS_SIM_SIM_DISCRETE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/sim/sim_common.h"

#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_target.h"

typedef struct
{
    int nx;
    int nu;
    int nz;
} sim_discrete_dims;



typedef struct
{
    external_function_generic *disc_dyn_fun;
    external_function_generic *disc_dyn_fun_jac;
    external_function_generic *disc_dyn_fun_jac_hess;
} discrete_model;


// workspace
typedef struct
{  
    double *rhs_forw_in;  // pointer to input vecotr for integration
    double *out_forw_traj;// pointer to output trajectory
} sim_discrete_workspace;

// memory
typedef struct
{
	// memory
	double time_sim;
	double time_ad;
	double time_la;

	// workspace structs

} sim_discrete_memory;



// discrete dims
acados_size_t sim_discrete_dims_calculate_size();
void *sim_discrete_dims_assign(void *config_, void *raw_memory);

// get & set functions
void sim_discrete_dims_set(void *config_, void *dims_, const char *field, const int *value);
void sim_discrete_dims_get(void *config_, void *dims_, const char *field, int* value);

// model
acados_size_t sim_discrete_model_calculate_size(void *config, void *dims_);
void *sim_discrete_model_assign(void *config, void *dims_, void *raw_memory);
int sim_discrete_model_set(void *model_, const char *field, void *value);

// opts
acados_size_t sim_discrete_opts_calculate_size(void *config, void *dims);
void *sim_discrete_opts_assign(void *config, void *dims, void *raw_memory);
void sim_discrete_opts_initialize_default(void *config, void *dims, void *opts_);
void sim_discrete_opts_update(void *config_, void *dims, void *opts_);
void sim_discrete_opts_set(void *config_, void *opts_, const char *field, void *value);


// memory
acados_size_t sim_discrete_memory_calculate_size(void *config, void *dims, void *opts_);
void *sim_discrete_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
int sim_discrete_memory_set(void *config_, void *dims_, void *mem_, const char *field, void *value);

// workspace
acados_size_t sim_discrete_workspace_calculate_size(void *config, void *dims, void *opts_);

// interface
void sim_discrete_config_initialize_default(void *config_);

// integrator
int sim_discrete(void *config, sim_in *in, sim_out *out, void *opts, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_DISCRETE_H_
