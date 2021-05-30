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


// standard
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/math.h"

#include "acados/sim/sim_collocation_utils.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_discrete.h"

// blasfeo
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_target.h"
// #include "blasfeo/include/blasfeo_d_aux_ext_dep.h" // can be included for printing while
// debugging



/************************************************
 * dims
 ************************************************/

acados_size_t sim_discrete_dims_calculate_size()
{
    acados_size_t size = sizeof(sim_discrete_dims);
    return size;
}

void *sim_discrete_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    sim_discrete_dims *dims = (sim_discrete_dims *) c_ptr;
    c_ptr += sizeof(sim_discrete_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;

    assert((char *) raw_memory + sim_discrete_dims_calculate_size() == c_ptr);
    return dims;
}



/************************************************
 * get & set functions
 ************************************************/

void sim_discrete_dims_set(void *config_, void *dims_, const char *field, const int *value)
{
    sim_discrete_dims *dims = dims_;

    if (!strcmp(field, "nx"))
    {
        dims->nx = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu = *value;
    }
    else if (!strcmp(field, "nz"))
    {
        dims->nz = *value;
    }
   else
    {
        printf("\nerror: sim_discrete_dims_set: field not available: %s\n", field);
        exit(1);
    }
}



void sim_discrete_dims_get(void *config_, void *dims_, const char *field, int *value)
{
    sim_discrete_dims *dims = dims_;

    if (!strcmp(field, "nx"))
    {
        *value = dims->nx;
    }
    else if (!strcmp(field, "nu"))
    {
        *value = dims->nu;
    }
    else if (!strcmp(field, "nz"))
    {
        *value = dims->nz;
    }
    else
    {
        printf("\nerror: sim_discrete_dims_get: field not available: %s\n", field);
        exit(1);
    }
}


/************************************************
 * opts
 ************************************************/

acados_size_t sim_discrete_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = NS_MAX;

    acados_size_t size = 0;

    size += sizeof(sim_opts);

    size += ns_max * ns_max * sizeof(double);  // A_mat
    size += ns_max * sizeof(double);           // b_vec
    size += ns_max * sizeof(double);           // c_vec

    size_t tmp0 = gauss_nodes_work_calculate_size(ns_max);
    size_t tmp1 = butcher_table_work_calculate_size(ns_max);
    acados_size_t work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    size += work_size;  // work

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}

void *sim_discrete_opts_assign(void *config_, void *dims, void *raw_memory)
{
    int ns_max = NS_MAX;

    char *c_ptr = (char *) raw_memory;

    sim_opts *opts = (sim_opts *) c_ptr;
    c_ptr += sizeof(sim_opts);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(ns_max * ns_max, &opts->A_mat, &c_ptr);
    assign_and_advance_double(ns_max, &opts->b_vec, &c_ptr);
    assign_and_advance_double(ns_max, &opts->c_vec, &c_ptr);

    // work
    acados_size_t tmp0 = gauss_nodes_work_calculate_size(ns_max);
    acados_size_t tmp1 = butcher_table_work_calculate_size(ns_max);
    acados_size_t work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    opts->work = c_ptr;
    c_ptr += work_size;

    assert((char *) raw_memory + sim_discrete_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}

void sim_discrete_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_discrete_dims *dims = (sim_discrete_dims *) dims_;
    sim_opts *opts = opts_;

    opts->ns = 3;  // GL 3
    int ns = opts->ns;

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    // gauss collocation nodes
    gauss_nodes(ns, opts->c_vec, opts->work);

    // butcher tableau
    butcher_table(ns, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    // default options
    opts->newton_iter = 3;
    opts->scheme = NULL;
    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
    opts->jac_reuse = true;
    opts->exact_z_output = false;

    // TODO(oj): check if constr h or cost depend on z, turn on in this case only.
    if (dims->nz > 0)
    {
        opts->output_z = true;
        opts->sens_algebraic = true;
    }
    else
    {
        opts->output_z = false;
        opts->sens_algebraic = false;
    }

    return;
}



void sim_discrete_opts_update(void *config_, void *dims, void *opts_)
{
    sim_opts *opts = opts_;

    int ns = opts->ns;

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    // gauss collocation nodes
    gauss_nodes(ns, opts->c_vec, opts->work);

    // butcher tableau
    butcher_table(ns, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    return;
}



void sim_discrete_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_set_(opts, field, value);
}



void sim_discrete_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_get_(config_, opts, field, value);
}



/************************************************
 * model
 ************************************************/

acados_size_t sim_discrete_model_calculate_size(void *config, void *dims_)
{
    acados_size_t size = 0;
    size += sizeof(discrete_model);

    return size;
}



void *sim_discrete_model_assign(void *config, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    discrete_model *data = (discrete_model *) c_ptr;
    c_ptr += sizeof(discrete_model);

    assert((char *) raw_memory + sim_discrete_model_calculate_size(config, dims) >= c_ptr);

    return data;
}



int sim_discrete_model_set(void *model_, const char *field, void *value)
{
    discrete_model *model = model_;

    if (!strcmp(field, "dyn_fun") || !strcmp(field, "disc_dyn_fun"))
    {
        model->disc_dyn_fun = value;
    }
    else if (!strcmp(field, "dyn_fun_jac") || !strcmp(field, "disc_dyn_fun_jac"))
    {
        model->disc_dyn_fun_jac = value;
    }
    else if (!strcmp(field, "dyn_fun_jac_hess") || !strcmp(field, "disc_dyn_fun_jac_hess"))
    {
        model->disc_dyn_fun_jac_hess = value;
    }
    else
    {
        printf("\nerror: sim_discrete_model_set: wrong field: %s\n", field);
        exit(1);
    }

    return ACADOS_SUCCESS;
}


/************************************************
 * memory
 ************************************************/

acados_size_t sim_discrete_memory_calculate_size(void *config, void *dims_, void *opts_)
{
    acados_size_t size = 0;

    size += sizeof(sim_discrete_memory);

    return size;
}



void *sim_discrete_memory_assign(void *config, void *dims_, void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_discrete_memory *mem = (sim_discrete_memory *) c_ptr;
    c_ptr += sizeof(sim_discrete_memory);

    return mem;
}



int sim_discrete_memory_set(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    printf("sim_discrete_memory_set field %s is not supported! \n", field);
    exit(1);
}



int sim_discrete_memory_set_to_zero(void *config_, void * dims_, void *opts_, void *mem_, const char *field)
{
    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "guesses"))
    {
        // no guesses/initialization in ERK
    }
    else
    {
        printf("sim_discrete_memory_set_to_zero field %s is not supported! \n", field);
        exit(1);
    }
}



void sim_discrete_memory_get(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    sim_discrete_memory *mem = mem_;

    if (!strcmp(field, "time_sim"))
    {
		double *ptr = value;
		*ptr = mem->time_sim;
	}
    else if (!strcmp(field, "time_sim_ad"))
    {
		double *ptr = value;
		*ptr = mem->time_ad;
	}
    else if (!strcmp(field, "time_sim_la"))
    {
		double *ptr = value;
		*ptr = mem->time_la;
	}
	else
	{
		printf("sim_discrete_memory_get field %s is not supported! \n", field);
		exit(1);
	}
}



/************************************************
 * workspace
 ************************************************/

acados_size_t sim_discrete_workspace_calculate_size(void *config, void *dims_, void *opts_)
{
    // ocp_nlp_dynamics_config *config = config_;
    sim_discrete_dims *dims = dims_;
    sim_opts *opts = opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    // int nx1 = dims->nx1;

    acados_size_t size = 0;

    size += sizeof(sim_discrete_workspace);

    return size;
}


/************************************************
 * functions
 ************************************************/

int sim_discrete_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    return ACADOS_SUCCESS;
}

static void *sim_discrete_cast_workspace(void *config, void *dims_, void *opts_, void *raw_memory)
{
    // typecast
    sim_discrete_dims *dims = (sim_discrete_dims *) dims_;
    sim_opts *opts = opts_;

    // dimension ints
    int nx = dims->nx;
    int nu = dims->nu;

    char *c_ptr = (char *) raw_memory;
    sim_discrete_workspace *workspace = (sim_discrete_workspace *) c_ptr;
    c_ptr += sizeof(sim_discrete_workspace);
    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + sim_discrete_workspace_calculate_size(config, dims_, opts) >= c_ptr);

    return (void *) workspace;
}



int sim_discrete(void *config, sim_in *in, sim_out *out, void *args, void *mem_, void *work_)
{
    acados_timer tot_timer, casadi_timer, la_timer;
    acados_tic(&tot_timer);

    // typecast
    sim_discrete_memory *mem = (sim_discrete_memory *) mem_;
    sim_opts *opts = (sim_opts *) args;
    sim_discrete_dims *dims = (sim_discrete_dims *) in->dims;
    discrete_model *model = in->model;
    sim_discrete_workspace *workspace =
        (sim_discrete_workspace *) sim_discrete_cast_workspace(config, dims, opts, work_);

    // necessary integers
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nz      = dims->nz;

    out->info->CPUtime = acados_toc(&tot_timer);

	mem->time_sim = out->info->CPUtime;
	mem->time_ad = out->info->ADtime;
	mem->time_la = out->info->LAtime;

    return ACADOS_SUCCESS;
}



void sim_discrete_config_initialize_default(void *config_)
{
    sim_config *config = config_;
    config->evaluate = &sim_discrete;
    config->precompute = &sim_discrete_precompute;
    // opts
    config->opts_calculate_size = &sim_discrete_opts_calculate_size;
    config->opts_assign = &sim_discrete_opts_assign;
    config->opts_initialize_default = &sim_discrete_opts_initialize_default;
    config->opts_update = &sim_discrete_opts_update;
    config->opts_set = &sim_discrete_opts_set;
    config->opts_get = &sim_discrete_opts_get;
    // memory & workspace
    config->memory_calculate_size = &sim_discrete_memory_calculate_size;
    config->memory_assign = &sim_discrete_memory_assign;
    config->memory_set = &sim_discrete_memory_set;
    config->memory_set_to_zero = &sim_discrete_memory_set_to_zero;
    config->memory_get = &sim_discrete_memory_get;
    config->workspace_calculate_size = &sim_discrete_workspace_calculate_size;
    // model
    config->model_calculate_size = &sim_discrete_model_calculate_size;
    config->model_assign = &sim_discrete_model_assign;
    config->model_set = &sim_discrete_model_set;
    // dims
    config->dims_calculate_size = &sim_discrete_dims_calculate_size;
    config->dims_assign = &sim_discrete_dims_assign;
    config->dims_set = &sim_discrete_dims_set;
    config->dims_get = &sim_discrete_dims_get;
    return;
}
