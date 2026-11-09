/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Copied from t8code example/advect/t8_advection.cxx at commit cda465320094e3179971ea5a09cfab6a35031407.
 * All functions and types defined here were renamed from t8_* to t8code_GPU_experiments_*. */

#include <sc_options.h>
#include <sc_statistics.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <advect/t8code_GPU_experiments_example_common.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_vector/t8_vec.hxx>
/* '#pragma omp requires unified_shared_memory' is in t8code_GPU_experiments_example_common.hxx */

#define MAX_FACES 8 /* The maximum number of faces of an element */
/* TODO: This is not memory efficient. If we run out of memory, we can optimize here. */
#define MAX_FACE_NEIGHBORS 4 /* The maximum number of face neighbors at one face in a 2:1 balanced forest */

/* Enum for statistics. */
typedef enum {
  ADVECT_ADAPT = 0,       /* adapt runtime */
  ADVECT_PARTITION,       /* partition runtime */
  ADVECT_PARTITION_PROCS, /* number of processes sent to in partition */
  ADVECT_PARTITION_DATA,  /* data partitioning runtime */
  ADVECT_BALANCE,         /* balance runtime */
  ADVECT_BALANCE_ROUNDS,  /* number of rounds in balance */
  ADVECT_GHOST,           /* ghost runtime */
  ADVECT_GHOST_SENT,      /* number of ghosts sent to other processes */
  ADVECT_GHOST_EXCHANGE,  /* ghost exchange runtime */
  ADVECT_GHOST_WAIT,      /* ghost exchange waittime */
  ADVECT_REPLACE,         /* forest_iterate_replace runtime */
  ADVECT_IO,              /* vtk runtime */
  ADVECT_ELEM_AVG,        /* average global number of elements (per time step) */
  ADVECT_INIT,            /* initialization runtime */
  ADVECT_AMR,             /* AMR runtime (adapt+partition+ghost+balance) including data exchange (partition/ghost) */
  ADVECT_NEIGHS,          /* neighbor finding runtime */
  ADVECT_FLUX,            /* flux computation runtime */
  ADVECT_DUMMY,           /* dummy operations to increase load (see -s option) */
  ADVECT_SOLVE,           /* solver runtime */
  ADVECT_TOTAL,           /* overall runtime */
  ADVECT_ERROR_INF,       /* l_infty error */
  ADVECT_ERROR_2,         /* L_2 error */
  ADVECT_VOL_LOSS,        /* The loss in volume (region with LS < 0) in percent */
  ADVECT_NUM_STATS        /* The number of statistics that we measure */
} advect_stats_t;

/* Names of statistics that we measure */
const char *advect_stat_names[ADVECT_NUM_STATS] = { "adapt",
                                                    "partition",
                                                    "partition_procs_sent",
                                                    "partition_data",
                                                    "balance",
                                                    "balance_rounds",
                                                    "ghost",
                                                    "ghost_sent",
                                                    "ghost_exchange",
                                                    "ghost_exchange_wait",
                                                    "replace",
                                                    "vtk_print",
                                                    "number_elements",
                                                    "init",
                                                    "AMR",
                                                    "neighbor_finding",
                                                    "flux_computation",
                                                    "dummy_ops",
                                                    "solve",
                                                    "total",
                                                    "l_infty_error",
                                                    "L_2",
                                                    "volume_loss_[%]" };

/** The per element data, stored as struct of arrays.
 *  Each array has one entry per local element (num_elements entries). */
typedef struct
{
  t8_locidx_t num_elements;                             /**< Number of local elements, length of each array */
  t8_3D_vec *midpoint;                                  /**< coordinates of element midpoint in R^3 */
  double *vol;                                          /**< Volume of the element */
  double *phi_new;                                      /**< Value of solution at midpoint in next time step */
  double (*fluxes)[MAX_FACES][MAX_FACE_NEIGHBORS];      /**< The fluxes to each neighbor at a given face.
                                                             Only the first SC_MAX (1, num_neighbors) entries are used. */
  int (*flux_valid)[MAX_FACES];                         /**< If > 0, this flux was computed, if <= 0 it was not computed.
                                                             -1 marks faces whose neighbors were not yet computed. */
  int *level;                                           /**< The refinement level of the element. */
  int *num_faces;                                       /**< The number of faces */
  int (*num_neighbors)[MAX_FACES];                      /**< Number of neighbors for each face */
  int (*dual_faces)[MAX_FACES][MAX_FACE_NEIGHBORS];     /**< The face indices of the neighbor elements.
                                                             Only the first num_neighbors entries are used. */
  t8_locidx_t (*neighs)[MAX_FACES][MAX_FACE_NEIGHBORS]; /**< Indices of the neighbor elements.
                                                             Only the first num_neighbors entries are used. */
  int8_t (*neigh_level)[MAX_FACES];                     /**< The level of the face neighbors at this face. */
  t8_3D_vec (*face_center)[MAX_FACES][MAX_FACE_NEIGHBORS]; /**< The center of each face. If a face has more than one
                                                                neighbor, one entry for each child face. Not used in 1D. */
  t8_3D_vec (*face_normal)[MAX_FACES][MAX_FACE_NEIGHBORS]; /**< The outward normal of each face (or child face). */
  double (*face_area)[MAX_FACES][MAX_FACE_NEIGHBORS];      /**< The area of each face (or child face). */
  double *dummy;                                           /**< Result of the dummy operations (see -s option).
                                                                Not used otherwise. */
} t8code_GPU_experiments_advect_element_data_t;

/* Allocate an array of n entries of type T.
 * We need this, since T8_ALLOC does not accept array types such as double[MAX_FACES]. */
template <typename T>
static void
t8code_GPU_experiments_advect_alloc_field (T *&field, const t8_locidx_t n)
{
  field = T8_ALLOC (T, n);
}

/* Allocate element data for num_elements elements. The entries are not initialized. */
static t8code_GPU_experiments_advect_element_data_t *
t8code_GPU_experiments_advect_element_data_new (const t8_locidx_t num_elements)
{
  t8code_GPU_experiments_advect_element_data_t *data = T8_ALLOC (t8code_GPU_experiments_advect_element_data_t, 1);

  data->num_elements = num_elements;
  t8code_GPU_experiments_advect_alloc_field (data->midpoint, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->vol, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->phi_new, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->fluxes, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->flux_valid, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->level, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->num_faces, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->num_neighbors, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->dual_faces, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->neighs, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->neigh_level, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->face_center, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->face_normal, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->face_area, num_elements);
  t8code_GPU_experiments_advect_alloc_field (data->dummy, num_elements);
  return data;
}

/* Free element data and set the pointer to NULL. */
static void
t8code_GPU_experiments_advect_element_data_destroy (t8code_GPU_experiments_advect_element_data_t **pdata)
{
  t8code_GPU_experiments_advect_element_data_t *data;

  T8_ASSERT (pdata != NULL);
  data = *pdata;
  T8_FREE (data->midpoint);
  T8_FREE (data->vol);
  T8_FREE (data->phi_new);
  T8_FREE (data->fluxes);
  T8_FREE (data->flux_valid);
  T8_FREE (data->level);
  T8_FREE (data->num_faces);
  T8_FREE (data->num_neighbors);
  T8_FREE (data->dual_faces);
  T8_FREE (data->neighs);
  T8_FREE (data->neigh_level);
  T8_FREE (data->face_center);
  T8_FREE (data->face_normal);
  T8_FREE (data->face_area);
  T8_FREE (data->dummy);
  T8_FREE (data);
  *pdata = NULL;
}

/* Copy all data of element isrc in src to element idst in dst.
 * The face geometry is not copied, since it is recomputed after each adapt and partition.
 * The dummy value is not copied, since it is never read. */
static void
t8code_GPU_experiments_advect_element_data_copy (t8code_GPU_experiments_advect_element_data_t *dst, const t8_locidx_t idst,
                                                 const t8code_GPU_experiments_advect_element_data_t *src,
                                                 const t8_locidx_t isrc)
{
  T8_ASSERT (0 <= idst && idst < dst->num_elements);
  T8_ASSERT (0 <= isrc && isrc < src->num_elements);
  dst->midpoint[idst] = src->midpoint[isrc];
  dst->vol[idst] = src->vol[isrc];
  dst->phi_new[idst] = src->phi_new[isrc];
  memcpy (dst->fluxes[idst], src->fluxes[isrc], sizeof (*dst->fluxes));
  memcpy (dst->flux_valid[idst], src->flux_valid[isrc], sizeof (*dst->flux_valid));
  dst->level[idst] = src->level[isrc];
  dst->num_faces[idst] = src->num_faces[isrc];
  memcpy (dst->num_neighbors[idst], src->num_neighbors[isrc], sizeof (*dst->num_neighbors));
  memcpy (dst->dual_faces[idst], src->dual_faces[isrc], sizeof (*dst->dual_faces));
  memcpy (dst->neighs[idst], src->neighs[isrc], sizeof (*dst->neighs));
  memcpy (dst->neigh_level[idst], src->neigh_level[isrc], sizeof (*dst->neigh_level));
}

/* Send the entries of one element data field from the partition of forest_from
 * to the partition of forest_to. Both arrays have one entry per local element. */
template <typename T>
static void
t8code_GPU_experiments_advect_partition_field (const t8_forest_t forest_from, const t8_forest_t forest_to, T *field_from,
                                               T *field_to)
{
  sc_array_t view_from, view_to;

  sc_array_init_data (&view_from, field_from, sizeof (T), t8_forest_get_local_num_leaf_elements (forest_from));
  sc_array_init_data (&view_to, field_to, sizeof (T), t8_forest_get_local_num_leaf_elements (forest_to));
  t8_forest_partition_data (forest_from, forest_to, &view_from, &view_to);
}

/** The description of the problem configuration.
 *  We store all necessary parameters, such as the initial level-set function, the flow function,
 *  data needed for adaptation etc.
 *  We also store the current forest and element data here.
 */
typedef struct
{
  t8code_GPU_experiments_flow_function_3d_fn u;      /**< Fluid field */
  int flow;                      /**< Number of the fluid field u, used to evaluate u on the device (see flow_eval) */
  t8code_GPU_experiments_example_level_set_fn phi_0; /**< Initial condition for phi */
  void *udata_for_phi;           /**< User data passed to phi */
  t8_forest_t forest;            /**< The forest in use */
  t8_forest_t forest_adapt;      /**< The forest after adaptation */
  t8code_GPU_experiments_advect_element_data_t *element_data; /**< The data of the local elements */
  t8code_GPU_experiments_advect_element_data_t
    *element_data_adapt; /**< element_data for the adapted forest, used during adaptation to interpolate values */
  /* We store the phi values in an extra array, since this data must exist for ghost
   * element as well and is communicated with other processes in ghost_exchange. */
  sc_array_t *phi_values;       /**< For each element and ghost its phi value. */
  sc_array_t *phi_values_adapt; /**< phi values for the adapted forest, used during adaptation to interpolate values. */
  sc_MPI_Comm comm;             /**< MPI communicator used */
  sc_statinfo_t stats[ADVECT_NUM_STATS]; /**< Runtimes and other statistics. */
  double t;                              /**< Current simulation time */
  double T;                              /**< End time */
  double cfl;                            /**< CFL number */
  double delta_t;                        /**< Current time step */
  double min_grad, max_grad;             /**< bounds for refinement */
  double min_vol;                        /**< minimum element volume at level 'level' */
  double band_width;                     /**< width of the refinement band */
  int num_time_steps;                    /**< Number of time steps computed so far.
                                        (If delta_t is constant then t = num_time_steps * delta_t) */
  int vtk_count;                         /**< If vtk output is enabled, count the number of pvtu files written. */
  int level;                             /**< Initial refinement level */
  int maxlevel;                          /**< Maximum refinement level */
  int volume_refine;                     /**< If >= refine elements only if their volume is greater
                                       than the minimum volume at level 'level + volume_refine' */
  int dim;                               /**< The dimension of the mesh */
  int dummy_op;                          /**< If true, we carry out more (but useless) operations
                                     per element, in order to simulate more computation load */
} t8code_GPU_experiments_advect_problem_t;

/* Return the phi value of a given local or ghost element.
 * 0 <= ielement < num_elements + num_ghosts
 */
static double
t8code_GPU_experiments_advect_element_get_phi (const t8code_GPU_experiments_advect_problem_t *problem, t8_locidx_t ielement)
{
  return *((double *) t8_sc_array_index_locidx (problem->phi_values, ielement));
}

/* Set the phi value of an element to a given entry */
static void
t8code_GPU_experiments_advect_element_set_phi (const t8code_GPU_experiments_advect_problem_t *problem, t8_locidx_t ielement, double phi)
{
  *((double *) t8_sc_array_index_locidx (problem->phi_values, ielement)) = phi;
}

/* Set the phi value of an element in the adapted forest to a given entry */
static void
t8code_GPU_experiments_advect_element_set_phi_adapt (const t8code_GPU_experiments_advect_problem_t *problem, t8_locidx_t ielement, double phi)
{
  *((double *) t8_sc_array_index_locidx (problem->phi_values_adapt, ielement)) = phi;
}

/* Adapt the forest. We refine if the level-set function is close to zero
 * and coarsen if it is larger than a given threshold. */
static int
t8code_GPU_experiments_advect_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id, const t8_eclass_t tree_class,
                 t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family,
                 [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  t8code_GPU_experiments_advect_problem_t *problem;
  double band_width, elem_diam;
  int level;
  t8_locidx_t offset;
  double phi;
  double vol_thresh;
  static int seed = 10000;

  srand (seed++);
  /* Get a pointer to the problem from the user data pointer of forest */
  problem = (t8code_GPU_experiments_advect_problem_t *) t8_forest_get_user_data (forest);
  /* Get the element's level */
  level = scheme->element_get_level (tree_class, elements[0]);
  if (level == problem->maxlevel && !is_family) {
    /* It is not possible to refine this level */
    return 0;
  }
  /* Compute the volume threshold. Elements larger than this and
   * close to the 0 level-set are refined */
  if (problem->volume_refine >= 0) {
    vol_thresh = problem->min_vol / (1 << (problem->dim * problem->volume_refine));
  }
  else {
    vol_thresh = 0;
  }
  /* Get the value of phi at this element */
  offset = t8_forest_get_tree_element_offset (forest_from, ltree_id);
  phi = t8code_GPU_experiments_advect_element_get_phi (problem, lelement_id + offset);

  /* Refine if close to levelset, coarsen if not */
  band_width = problem->band_width;
  elem_diam = t8_forest_element_diam (forest_from, ltree_id, elements[0]);
  if (fabs (phi) > 2 * band_width * elem_diam) {
    /* coarsen if this is a family and level is not too small */
    return -(is_family && level > problem->level);
  }
  else if (fabs (phi) < band_width * elem_diam && problem->element_data->vol[lelement_id + offset] > vol_thresh) {
    /* refine if level is not too large */
    return level < problem->maxlevel;
  }
  return 0;
}

/* Compute the total volume of the elements with negative phi value */
static double
t8code_GPU_experiments_advect_level_set_volume (const t8code_GPU_experiments_advect_problem_t *problem)
{
  t8_locidx_t num_local_elements, ielem;
  const t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;
  double volume = 0, global_volume = 0;
  double phi;

  num_local_elements = t8_forest_get_local_num_leaf_elements (problem->forest);

  for (ielem = 0; ielem < num_local_elements; ielem++) {
    phi = t8code_GPU_experiments_advect_element_get_phi (problem, ielem);
    if (phi < 0) {
      volume += data->vol[ielem];
    }
  }
  sc_MPI_Allreduce (&volume, &global_volume, 1, sc_MPI_DOUBLE, sc_MPI_SUM, problem->comm);
  return global_volume;
}

/* Compute the relative l_infty error of the stored phi values compared to a
 * given analytical function at time problem->t */
static double
t8code_GPU_experiments_advect_l_infty_rel (const t8code_GPU_experiments_advect_problem_t *problem, t8code_GPU_experiments_example_level_set_fn analytical_sol, double distance)
{
  t8_locidx_t num_local_elements, ielem;
  const t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;
  double phi;
  double ana_sol;
  double error[2] = { -1, 0 }, el_error, global_error[2];

  num_local_elements = t8_forest_get_local_num_leaf_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    /* Compute the analytical solution */
    ana_sol = analytical_sol (data->midpoint[ielem], problem->t, problem->udata_for_phi);
    if (fabs (ana_sol) < distance) {
      /* Compute the error as the stored value at the midpoint of this element
       * minus the solution at this midpoint */
      phi = t8code_GPU_experiments_advect_element_get_phi (problem, ielem);
      el_error = fabs ((phi - ana_sol));
      error[0] = SC_MAX (error[0], el_error);
      /* Compute the l_infty norm of the analytical solution */
      error[1] = SC_MAX (error[1], ana_sol);
    }
  }
  /* Compute the maximum of the error among all processes */
  sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_MAX, problem->comm);

  /* Return the relative error, that is the l_infty error divided by
   * the l_infty norm of the analytical solution */
  return global_error[0] / global_error[1];
}

static double
t8code_GPU_experiments_advect_l_2_rel (const t8code_GPU_experiments_advect_problem_t *problem, t8code_GPU_experiments_example_level_set_fn analytical_sol, double distance)
{
  t8_locidx_t num_local_elements, ielem, count = 0;
  const t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;
  double phi;
  double diff, ana_sol;
  double error[2] = { 0, 0 }, el_error, global_error[2];

  num_local_elements = t8_forest_get_local_num_leaf_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    /* Compute the analytical solution */
    ana_sol = analytical_sol (data->midpoint[ielem], problem->t, problem->udata_for_phi);
    if (fabs (ana_sol) < distance) {
      count++;
      /* Compute the error as the stored value at the midpoint of this element minus the solution at this midpoint */
      phi = t8code_GPU_experiments_advect_element_get_phi (problem, ielem);
      diff = fabs (phi - ana_sol);
      el_error = diff * diff * data->vol[ielem];
      error[0] += el_error;
      error[1] += ana_sol * ana_sol * data->vol[ielem];
    }
  }
  t8_debugf ("[advect] L_2 %e  %e\n", error[0], error[1]);
  t8_debugf ("[advect] L_2 %i elems\n", count);
  /* Compute the maximum of the error among all processes */
  sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_SUM, problem->comm);

  /* Return the relative error, that is the l_infty error divided by the l_infty norm of the analytical solution */
  return sqrt (global_error[0]) / sqrt (global_error[1]);
}

/* The following functions are called inside the element loop of the solver, which can run on the device.
 * They must not call t8code or sc functions and only work on precomputed data. */
#pragma omp declare target

/* Evaluate the flow function with number flow (see t8code_GPU_experiments_advect_choose_flow)
 * at the point x and time t and store the result in u.
 * We cannot call problem->u on the device, since it is a host function pointer. */
static void
t8code_GPU_experiments_advect_flow_eval (const int flow, const t8_3D_vec &x, const double t, t8_3D_vec &u)
{
  switch (flow) {
  case 1:
    t8code_GPU_experiments_flow_constant_one_x_vec (x, t, u);
    break;
  case 2:
    t8code_GPU_experiments_flow_constant_one_xyz_vec (x, t, u);
    break;
  case 3:
    t8code_GPU_experiments_flow_incomp_cube_flow (x, t, u);
    break;
  case 4:
    t8code_GPU_experiments_flow_rotation_2d (x, t, u);
    break;
  case 5:
    t8code_GPU_experiments_flow_around_circle (x, t, u);
    break;
  case 6:
    t8code_GPU_experiments_flow_stokes_flow_sphere_shell (x, t, u);
    break;
  case 7:
    t8code_GPU_experiments_flow_around_circle_with_angular_velocity (x, t, u);
    break;
  default:
    /* Not reachable, flow is checked in main */
    u[0] = u[1] = u[2] = 0;
  }
}

/* Compute the flux across face number face of a 1D element with given midpoint and vol.
 * phi_plus is the phi value of the element, phi_minus the phi value of the neighbor
 * (0 if the element is at the boundary). */
static double
t8code_GPU_experiments_advect_flux_upwind_1d (const int flow, const double t, const t8_3D_vec &midpoint, const double vol,
                                              const double phi_plus, const double phi_minus, const int face)
{
  t8_3D_vec x_j_half;
  int idim;
  t8_3D_vec u_at_x_j_half;
  double phi;
  int sign;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          x_j_half
   */
  /* Compute x_j_half */
  for (idim = 0; idim < 3; idim++) {
    x_j_half[idim] = (midpoint[idim] - (idim == 0 ? vol / 2 : 0));
  }
  /* Compute u at the interval boundary. */
  t8code_GPU_experiments_advect_flow_eval (flow, x_j_half, t, u_at_x_j_half);
  /* In 1D we are only interested in the first coordinate of u */

  sign = face == 0 ? -1 : 1;
  if (sign * u_at_x_j_half[0] >= 0) {
    /* we have outflow */
    phi = -phi_plus;
  }
  else {
    /* we have inflow */
    phi = phi_minus;
  }
  return u_at_x_j_half[0] * phi;
}

/* Compute the flux across a face with given center, outward normal and area.
 * phi_plus is the phi value of the element, phi_minus the phi value of the neighbor. */
static double
t8code_GPU_experiments_advect_flux_upwind (const int flow, const double t, const double phi_plus, const double phi_minus,
                                           const t8_3D_vec &face_center, const t8_3D_vec &normal, const double area)
{
  t8_3D_vec u_at_face_center;
  double normal_times_u;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          face_center
   */

  /* Compute u at the face center. */
  t8code_GPU_experiments_advect_flow_eval (flow, face_center, t, u_at_face_center);
  /* Compute the dot-product of u and the normal vector */
  normal_times_u = t8_dot (normal, u_at_face_center);

  if (normal_times_u >= 0) {
    return -phi_plus * normal_times_u * area;
  }
  else {
    /* u flows into the element */
    return -phi_minus * normal_times_u * area;
  }
}

/* Compute all fluxes of element lelement and store them in data->fluxes.
 * The element only writes its own fluxes. Thus, each inner face is computed from both sides,
 * but elements can be processed in parallel.
 * phi holds the phi values of all local elements and ghosts. */
static void
t8code_GPU_experiments_advect_compute_fluxes (t8code_GPU_experiments_advect_element_data_t *data,
                                              const t8_locidx_t lelement, const double *phi, const int dim,
                                              const int flow, const double t)
{
  const double phi_plus = phi[lelement];
  const int num_faces = data->num_faces[lelement];
  int iface, ineigh;

  for (iface = 0; iface < num_faces; iface++) {
    const int num_neighbors = data->num_neighbors[lelement][iface];
    if (dim == 1) {
      /* At the boundary we use phi = 0 */
      const double phi_minus = num_neighbors > 0 ? phi[data->neighs[lelement][iface][0]] : 0;
      data->fluxes[lelement][iface][0] = t8code_GPU_experiments_advect_flux_upwind_1d (
        flow, t, data->midpoint[lelement], data->vol[lelement], phi_plus, phi_minus, iface);
    }
    else if (num_neighbors <= 1) {
      /* There is exactly one face neighbor or the element is at the domain boundary.
       * At the boundary we enforce outflow boundary conditions by using the phi value of the element itself. */
      const double phi_minus = num_neighbors == 1 ? phi[data->neighs[lelement][iface][0]] : phi_plus;
      data->fluxes[lelement][iface][0]
        = t8code_GPU_experiments_advect_flux_upwind (flow, t, phi_plus, phi_minus, data->face_center[lelement][iface][0],
                                                     data->face_normal[lelement][iface][0],
                                                     data->face_area[lelement][iface][0]);
    }
    else {
      /* The face is hanging, we compute the flux to each neighbor across the corresponding child face. */
      for (ineigh = 0; ineigh < num_neighbors; ineigh++) {
        const double phi_minus = phi[data->neighs[lelement][iface][ineigh]];
        data->fluxes[lelement][iface][ineigh] = t8code_GPU_experiments_advect_flux_upwind (
          flow, t, phi_plus, phi_minus, data->face_center[lelement][iface][ineigh],
          data->face_normal[lelement][iface][ineigh], data->face_area[lelement][iface][ineigh]);
      }
    }
  }
}

/* Simulate more load per element (see -s option).
 * We store the result in data->dummy, which is not used otherwise. */
static void
t8code_GPU_experiments_advect_dummy_op (t8code_GPU_experiments_advect_element_data_t *data, const t8_locidx_t lelement)
{
  int i, j;
  double dummy = 0;

  for (i = 1; i < 5; i++) {
    dummy *= i;
    for (j = 0; j < 5; j++) {
      dummy += pow (i, j);
    }
  }
  data->dummy[lelement] = dummy;
}

/* Sum up the fluxes of element lelement and compute its phi value in the next time step. */
static void
t8code_GPU_experiments_advect_advance_element (t8code_GPU_experiments_advect_element_data_t *data,
                                               const t8_locidx_t lelement, const double *phi, const double delta_t)
{
  int iface, ineigh;
  double flux_sum = 0;
  int num_faces;

  num_faces = data->num_faces[lelement];
  /* Sum all the fluxes */
  for (iface = 0; iface < num_faces; iface++) {
    for (ineigh = 0; ineigh < SC_MAX (1, data->num_neighbors[lelement][iface]); ineigh++) {
      flux_sum += data->fluxes[lelement][iface][ineigh];
    }
  }
  /* Phi^t = dt/dx * (f_(j-1/2) - f_(j+1/2)) + Phi^(t-1) */
  data->phi_new[lelement] = (delta_t / data->vol[lelement]) * flux_sum + phi[lelement];
}

#pragma omp end declare target

/* Compute element midpoint and vol and store them for element ielem in data. */
static void
t8code_GPU_experiments_advect_compute_element_data (t8code_GPU_experiments_advect_problem_t *problem,
                                                    t8code_GPU_experiments_advect_element_data_t *data,
                                                    const t8_locidx_t ielem, const t8_element_t *element,
                                                    const t8_locidx_t ltreeid)
{
  /* Compute the midpoint coordinates of element */
  t8_forest_element_centroid (problem->forest, ltreeid, element, data->midpoint[ielem].data ());
  /* Compute the length of this element */
  data->vol[ielem] = t8_forest_element_volume (problem->forest, ltreeid, element);
}

/* Compute the face neighbors of element ielem across face iface and store their
 * indices, dual faces and level in data. */
static void
t8code_GPU_experiments_advect_compute_face_neighbors (const t8code_GPU_experiments_advect_problem_t *problem,
                                                      t8code_GPU_experiments_advect_element_data_t *data,
                                                      const t8_locidx_t ielem, const t8_element_t *element,
                                                      const t8_locidx_t ltreeid, const int iface)
{
  const t8_scheme *scheme = t8_forest_get_scheme (problem->forest);
  const t8_element_t **neighbors;
  int *dual_faces;
  t8_locidx_t *neigh_indices;
  t8_eclass_t neigh_eclass;
  int num_neighbors, ineigh;

  /* This allocates neighbors, dual_faces and neigh_indices if num_neighbors > 0 */
  t8_forest_leaf_face_neighbors (problem->forest, ltreeid, element, &neighbors, iface, &dual_faces, &num_neighbors,
                                 &neigh_indices, &neigh_eclass);
  SC_CHECK_ABORTF (num_neighbors <= MAX_FACE_NEIGHBORS, "Element has %i neighbors at face %i, at most %i are supported.\n",
                   num_neighbors, iface, MAX_FACE_NEIGHBORS);
  data->num_neighbors[ielem][iface] = num_neighbors;
  for (ineigh = 0; ineigh < num_neighbors; ineigh++) {
    data->neighs[ielem][iface][ineigh] = neigh_indices[ineigh];
    data->dual_faces[ielem][iface][ineigh] = dual_faces[ineigh];
    data->neigh_level[ielem][iface] = scheme->element_get_level (neigh_eclass, neighbors[ineigh]);
  }
  if (num_neighbors > 0) {
    T8_FREE (neighbors);
    T8_FREE (dual_faces);
    T8_FREE (neigh_indices);
  }
}

/* Compute the center, outward normal and area of face iface of element ielem and store them in data.
 * If the face has more than one neighbor (hanging face), we store these values for each child face
 * at the face instead, in the same order as the neighbors.
 * The face neighbors of ielem must be computed before. */
static void
t8code_GPU_experiments_advect_compute_face_geometry (const t8code_GPU_experiments_advect_problem_t *problem,
                                                     t8code_GPU_experiments_advect_element_data_t *data,
                                                     const t8_locidx_t ielem, const t8_element_t *element,
                                                     const t8_locidx_t ltreeid, const int iface)
{
  const int num_neighbors = data->num_neighbors[ielem][iface];

  if (num_neighbors <= 1) {
    t8_forest_element_face_centroid (problem->forest, ltreeid, element, iface,
                                     data->face_center[ielem][iface][0].data ());
    t8_forest_element_face_normal (problem->forest, ltreeid, element, iface, data->face_normal[ielem][iface][0].data ());
    data->face_area[ielem][iface][0] = t8_forest_element_face_area (problem->forest, ltreeid, element, iface);
  }
  else {
    const t8_scheme *scheme = t8_forest_get_scheme (problem->forest);
    const t8_eclass_t eclass = t8_forest_get_tree_class (problem->forest, ltreeid);
    t8_element_t **face_children;
    int ichild, child_face;

    /* Compute the children of the element at the face */
    const int num_face_children = scheme->element_get_num_face_children (eclass, element, iface);
    T8_ASSERT (num_face_children == num_neighbors);
    face_children = T8_ALLOC (t8_element_t *, num_face_children);
    scheme->element_new (eclass, num_face_children, face_children);
    scheme->element_get_children_at_face (eclass, element, iface, face_children, num_face_children, NULL);
    for (ichild = 0; ichild < num_face_children; ichild++) {
      child_face = scheme->element_face_get_child_face (eclass, element, iface, ichild);
      t8_forest_element_face_centroid (problem->forest, ltreeid, face_children[ichild], child_face,
                                       data->face_center[ielem][iface][ichild].data ());
      t8_forest_element_face_normal (problem->forest, ltreeid, face_children[ichild], child_face,
                                     data->face_normal[ielem][iface][ichild].data ());
      data->face_area[ielem][iface][ichild]
        = t8_forest_element_face_area (problem->forest, ltreeid, face_children[ichild], child_face);
    }
    /* clean-up */
    scheme->element_destroy (eclass, num_face_children, face_children);
    T8_FREE (face_children);
  }
}

/* Compute the face neighbors and, if dim > 1, the face geometry of all local elements.
 * This calls t8code and must be done on the host each time the forest changed. */
static void
t8code_GPU_experiments_advect_compute_mesh_data (const t8code_GPU_experiments_advect_problem_t *problem)
{
  t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;
  t8_locidx_t itree, ielement, lelement;
  int iface;

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, lelement = 0; itree < num_local_trees; itree++) {
    const t8_locidx_t num_tree_elements = t8_forest_get_tree_num_leaf_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_tree_elements; ielement++, lelement++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (problem->forest, itree, ielement);
      for (iface = 0; iface < data->num_faces[lelement]; iface++) {
        t8code_GPU_experiments_advect_compute_face_neighbors (problem, data, lelement, element, itree, iface);
        if (problem->dim > 1) {
          t8code_GPU_experiments_advect_compute_face_geometry (problem, data, lelement, element, itree, iface);
        }
      }
    }
  }
}

/* Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the phi value of its parent.
 * If elements are coarsened, the parent gets the average phi value of the children.
 */
/* outgoing are the old elements and incoming the new ones */
/* TODO: If coarsening, weight the phi values by volume of the children:
 *       phi_E = sum (phi_Ei *vol(E_i)/vol(E))
 *       Similar formula for refining?
 */
static void
t8code_GPU_experiments_advect_replace (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                   const t8_scheme *scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                   t8_locidx_t first_incoming)
{
  t8code_GPU_experiments_advect_problem_t *problem;
  t8code_GPU_experiments_advect_element_data_t *data_in, *data_out;
  t8_locidx_t first_incoming_data, first_outgoing_data;
  int i, iface;
  double phi_old;

  /* Get the problem description */

  problem = (t8code_GPU_experiments_advect_problem_t *) t8_forest_get_user_data (forest_new);
  T8_ASSERT (forest_old == problem->forest);
  T8_ASSERT (forest_new == problem->forest_adapt);
  /* Get the element data and the indices of the elements in it */
  first_incoming_data = first_incoming + t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing_data = first_outgoing + t8_forest_get_tree_element_offset (forest_old, which_tree);
  data_out = problem->element_data;
  data_in = problem->element_data_adapt;

  /* Get the old phi value (used in the cases with num_outgoing = 1) */
  phi_old = t8code_GPU_experiments_advect_element_get_phi (problem, first_outgoing_data);
  if (refine == 0) {
    T8_ASSERT (num_incoming == num_outgoing && num_incoming == 1);
    /* The element is not changed, copy phi and vol */
    t8code_GPU_experiments_advect_element_data_copy (data_in, first_incoming_data, data_out, first_outgoing_data);
    t8code_GPU_experiments_advect_element_set_phi_adapt (problem, first_incoming_data, phi_old);
#if T8_ENABLE_DEBUG
    /* Get a pointer to the new element */
    const t8_element_t *element
      = t8_forest_get_leaf_element_in_tree (problem->forest_adapt, which_tree, first_incoming);
    /* Debug check number of faces */
    T8_ASSERT (data_in->num_faces[first_incoming_data] == scheme->element_get_num_faces (tree_class, element));
#endif
    /* Set the neighbor entries to uninitialized */
    for (iface = 0; iface < data_in->num_faces[first_incoming_data]; iface++) {
      data_in->num_neighbors[first_incoming_data][iface] = 0;
      data_in->flux_valid[first_incoming_data][iface] = -1;
    }
  }
  else if (refine == 1) {
    T8_ASSERT (num_outgoing == 1);
#if T8_ENABLE_DEBUG
    /* Ensure that the number of incoming elements matches the
     * number of children of the outgoing element. */
    const t8_element_t *element_outgoing = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing);
    const int num_children = scheme->element_get_num_children (tree_class, element_outgoing);
    T8_ASSERT (num_incoming == num_children);
#endif
    /* The old element is refined, we copy the phi values and compute the new midpoints */
    for (i = 0; i < num_incoming; i++) {
      /* Get a pointer to the new element */
      const t8_element_t *element
        = t8_forest_get_leaf_element_in_tree (problem->forest_adapt, which_tree, first_incoming + i);
      /* Compute midpoint and vol of the new element */
      t8code_GPU_experiments_advect_compute_element_data (problem, data_in, first_incoming_data + i, element, which_tree);
      t8code_GPU_experiments_advect_element_set_phi_adapt (problem, first_incoming_data + i, phi_old);
      /* Set the neighbor entries to uninitialized */
      const int num_new_faces = scheme->element_get_num_faces (tree_class, element);
      data_in->num_faces[first_incoming_data + i] = num_new_faces;
      for (iface = 0; iface < num_new_faces; iface++) {
        data_in->num_neighbors[first_incoming_data + i][iface] = 0;
        data_in->flux_valid[first_incoming_data + i][iface] = -1;
      }
      /* Update the level */
      data_in->level[first_incoming_data + i] = data_out->level[first_outgoing_data] + 1;
    }
  }
  else {
    double phi = 0;
    T8_ASSERT (refine = -1);
    T8_ASSERT (num_incoming == 1);
#if T8_ENABLE_DEBUG
    /* Ensure that the number of outgoing elements matches the
     * number of siblings of the first outgoing element. */
    const t8_element_t *element_outgoing = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing);
    const int num_siblings = scheme->element_get_num_siblings (tree_class, element_outgoing);
    T8_ASSERT (num_outgoing == num_siblings);
#endif
    /* The old elements form a family which is coarsened. We compute the average
     * phi value and set it as the new phi value */
    /* Get a pointer to the new element */
    const t8_element_t *element
      = t8_forest_get_leaf_element_in_tree (problem->forest_adapt, which_tree, first_incoming);
    /* Compute midpoint and vol of the new element */
    t8code_GPU_experiments_advect_compute_element_data (problem, data_in, first_incoming_data, element, which_tree);

    /* Compute average of phi */
    for (i = 0; i < num_outgoing; i++) {
      phi += t8code_GPU_experiments_advect_element_get_phi (problem, first_outgoing_data + i);
    }
    phi /= num_outgoing;
    t8code_GPU_experiments_advect_element_set_phi_adapt (problem, first_incoming_data, phi);
    /* Set the neighbor entries to uninitialized */
    data_in->num_faces[first_incoming_data] = data_out->num_faces[first_outgoing_data];
    T8_ASSERT (data_in->num_faces[first_incoming_data] == scheme->element_get_num_faces (tree_class, element));
    for (iface = 0; iface < data_in->num_faces[first_incoming_data]; iface++) {
      data_in->num_neighbors[first_incoming_data][iface] = 0;
      data_in->flux_valid[first_incoming_data][iface] = -1;
    }
    /* update the level */
    data_in->level[first_incoming_data] = data_out->level[first_outgoing_data] - 1;
  }
}

/* Invalidate the neighbor and flux information of all local elements */
static void
t8code_GPU_experiments_advect_problem_elements_invalidate (t8code_GPU_experiments_advect_problem_t *problem)
{

  t8_locidx_t lelement, num_local_elem;
  int iface;
  t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;

  num_local_elem = t8_forest_get_local_num_leaf_elements (problem->forest);
  T8_ASSERT (num_local_elem <= data->num_elements);
  /* invalidate all elements */
  for (lelement = 0; lelement < num_local_elem; lelement++) {
    for (iface = 0; iface < data->num_faces[lelement]; iface++) {
      /* TODO: make face number dim independent */
      if (data->num_neighbors[lelement][iface] > 0) {
        data->num_neighbors[lelement][iface] = 0;
        data->flux_valid[lelement][iface] = -1;
      }
    }
  }
}

/* Adapt the forest and interpolate the phi values to the new grid,
 * compute the new u values on the grid */
static void
t8code_GPU_experiments_advect_problem_adapt (t8code_GPU_experiments_advect_problem_t *problem, int measure_time)
{
  t8_locidx_t num_elems_p_ghosts, num_elems;
  double adapt_time, balance_time = 0, ghost_time;
  t8_locidx_t ghost_sent;
  int balance_rounds, did_balance = 0;
  double replace_time;

  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&problem->forest_adapt);
  /* Enable profiling to measure the runtime */
  t8_forest_set_profiling (problem->forest_adapt, 1);
  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (problem->forest_adapt, problem);
  /* Set the adapt function */
  t8_forest_set_adapt (problem->forest_adapt, problem->forest, t8code_GPU_experiments_advect_adapt, 0);
  if (problem->maxlevel - problem->level > 1) {
    /* We also want to balance the forest if the difference in refinement levels is greater 1 */
    t8_forest_set_balance (problem->forest_adapt, NULL, 1);
    did_balance = 1;
  }
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost (problem->forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (problem->forest_adapt);

  /* Store the runtimes in problems stats */
  if (measure_time) {
    adapt_time = t8_forest_profile_get_adapt_time (problem->forest_adapt);
    ghost_time = t8_forest_profile_get_ghost_time (problem->forest_adapt, &ghost_sent);
    if (did_balance) {
      balance_time = t8_forest_profile_get_balance_time (problem->forest_adapt, &balance_rounds);
    }
    sc_stats_accumulate (&problem->stats[ADVECT_ADAPT], adapt_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST], ghost_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_SENT], ghost_sent);
    if (did_balance) {
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE], balance_time);
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE_ROUNDS], balance_rounds);
    }
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_ADAPT].count = 1;
    problem->stats[ADVECT_BALANCE].count = 1;
    problem->stats[ADVECT_GHOST].count = 1;
    problem->stats[ADVECT_GHOST_SENT].count = 1;
  }

  /* Allocate new memory for the element_data of the advected forest */
  num_elems = t8_forest_get_local_num_leaf_elements (problem->forest_adapt);
  num_elems_p_ghosts = num_elems + t8_forest_get_num_ghosts (problem->forest_adapt);
  problem->element_data_adapt = t8code_GPU_experiments_advect_element_data_new (num_elems);
  problem->phi_values_adapt = sc_array_new_count (sizeof (double), num_elems_p_ghosts);
  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without repartitioning. */
  replace_time = -sc_MPI_Wtime ();
  t8_forest_iterate_replace (problem->forest_adapt, problem->forest, t8code_GPU_experiments_advect_replace);
  replace_time += sc_MPI_Wtime ();
  if (measure_time) {
    sc_stats_accumulate (&problem->stats[ADVECT_REPLACE], replace_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR], ghost_time + adapt_time + balance_time + replace_time);
    problem->stats[ADVECT_REPLACE].count = 1;
    problem->stats[ADVECT_AMR].count = 1;
  }
  /* clean the old element data */
  t8code_GPU_experiments_advect_problem_elements_invalidate (problem);
  t8code_GPU_experiments_advect_element_data_destroy (&problem->element_data);
  sc_array_destroy (problem->phi_values);
  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = problem->forest_adapt;
  problem->forest_adapt = NULL;
  /* Set the elem data to the adapted elem data */
  problem->element_data = problem->element_data_adapt;
  problem->element_data_adapt = NULL;
  /* Set the phi values to the adapted phi values */
  problem->phi_values = problem->phi_values_adapt;
  problem->phi_values_adapt = NULL;
}

/* Re-partition the forest and element data of a problem */
static void
t8code_GPU_experiments_advect_problem_partition (t8code_GPU_experiments_advect_problem_t *problem, int measure_time)
{
  t8_forest_t forest_partition;
  sc_array_t phi_view, phi_view_new;
  t8code_GPU_experiments_advect_element_data_t *data, *new_data;
  sc_array_t *new_phi;
  t8_locidx_t num_local_elements, num_local_elements_new;
  t8_locidx_t num_ghosts_new;
  int procs_sent;
  t8_locidx_t ghost_sent;
  double partition_time, ghost_time;

  /* Partition the forest and create its ghost layer */
  /* ref the current forest, since we still need access to it */
  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_partition);
  /* Enable profiling to measure runtime */
  t8_forest_set_profiling (forest_partition, 1);
  /* Partition the forest and create ghosts */
  t8_forest_set_partition (forest_partition, problem->forest, 0);
  t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_partition);

  /* Add runtimes to internal stats */
  if (measure_time) {
    partition_time = t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
    ghost_time = t8_forest_profile_get_ghost_time (forest_partition, &ghost_sent);
    sc_stats_accumulate (&problem->stats[ADVECT_PARTITION], partition_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST], ghost_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR], ghost_time + partition_time);
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_PARTITION].count = 1;
    problem->stats[ADVECT_GHOST].count = 1;
  }
  /* Partition the data */
  num_local_elements = t8_forest_get_local_num_leaf_elements (problem->forest);
  num_local_elements_new = t8_forest_get_local_num_leaf_elements (forest_partition);
  num_ghosts_new = t8_forest_get_num_ghosts (forest_partition);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&phi_view, problem->phi_values, 0, num_local_elements);
  /* Allocate the data arrays for the partitioned elements */
  data = problem->element_data;
  new_data = t8code_GPU_experiments_advect_element_data_new (num_local_elements_new);
  new_phi = sc_array_new_count (sizeof (double), num_local_elements_new + num_ghosts_new);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&phi_view_new, new_phi, 0, num_local_elements_new);
  /* Perform the data partition, one call per element data field */
  partition_time = -sc_MPI_Wtime ();
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->midpoint, new_data->midpoint);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->vol, new_data->vol);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->phi_new, new_data->phi_new);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->fluxes, new_data->fluxes);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->flux_valid,
                                                 new_data->flux_valid);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->level, new_data->level);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->num_faces, new_data->num_faces);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->num_neighbors,
                                                 new_data->num_neighbors);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->dual_faces,
                                                 new_data->dual_faces);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->neighs, new_data->neighs);
  t8code_GPU_experiments_advect_partition_field (problem->forest, forest_partition, data->neigh_level,
                                                 new_data->neigh_level);
  /* The face geometry is not partitioned, since it is recomputed after each partition.
   * The dummy value is not partitioned, since it is never read. */
  t8_forest_partition_data (problem->forest, forest_partition, &phi_view, &phi_view_new);
  partition_time += sc_MPI_Wtime ();
  if (measure_time) {
    sc_stats_accumulate (&problem->stats[ADVECT_PARTITION_DATA], partition_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR], partition_time);
    problem->stats[ADVECT_PARTITION_DATA].count = 1;
    problem->stats[ADVECT_AMR].count = 1;
    t8_debugf ("status ghost: %f\n\n", ghost_time);
  }

  /* destroy the old forest and the element data */
  t8code_GPU_experiments_advect_problem_elements_invalidate (problem);
  t8_forest_unref (&problem->forest);
  problem->forest = forest_partition;
  t8code_GPU_experiments_advect_element_data_destroy (&problem->element_data);
  problem->element_data = new_data;
  sc_array_destroy (problem->phi_values);
  problem->phi_values = new_phi;
}

static t8_cmesh_t
t8code_GPU_experiments_advect_create_cmesh (sc_MPI_Comm comm, int cube_type, const char *mshfile, int level, int dim, int use_cad_geometry)
{
  if (mshfile != NULL) {
    /* Load from .msh file and partition */
    t8_cmesh_t cmesh, cmesh_partition;
    T8_ASSERT (mshfile != NULL);

    t8_cmesh_init (&cmesh);
    t8_cmesh_from_msh_file (&cmesh, mshfile, 0, comm, dim, 0, use_cad_geometry);
    /* The partitioning of the cad geometry is not yet available */
    if (use_cad_geometry) {
      t8_productionf ("cmesh was not partitioned. Partitioning is not yet "
                      "available with the curved geometry\n");
      return cmesh;
    }
    /* partition this cmesh according to the initial refinement level */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_partition_uniform (cmesh_partition, level, t8_scheme_new_default ());
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_commit (cmesh_partition, comm);
    return cmesh_partition;
  }
  else {
    t8_cmesh_t cmesh;
    if (cube_type == 7) {
      t8_cmesh_init (&cmesh);
      t8_cmesh_new_periodic_hybrid (cmesh, comm);
      return cmesh;
    }
    else if (cube_type == 8) {
      t8_cmesh_init (&cmesh);
      t8_cmesh_new_hypercube_hybrid (cmesh, comm, 1);
      return cmesh;
    }
    else if (cube_type == 9) {
      t8_cmesh_init (&cmesh);
      t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_PYRAMID, comm, 0, 0, 0);
      return cmesh;
    }
    else {
      T8_ASSERT (T8_ECLASS_ZERO <= cube_type && cube_type < T8_ECLASS_COUNT);
      t8_cmesh_init (&cmesh);
      t8_cmesh_new_hypercube (&cmesh, (t8_eclass_t) cube_type, comm, 0, 0, 1);
      return cmesh;
    }
  }
}

static t8code_GPU_experiments_flow_function_3d_fn
t8code_GPU_experiments_advect_choose_flow (int flow_arg)
{
  switch (flow_arg) {
  case 1:
    return t8code_GPU_experiments_flow_constant_one_x_vec;
  case 2:
    return t8code_GPU_experiments_flow_constant_one_xyz_vec;
  case 3:
    return t8code_GPU_experiments_flow_incomp_cube_flow;
  case 4:
    return t8code_GPU_experiments_flow_rotation_2d;
  case 5:
    return t8code_GPU_experiments_flow_around_circle;
  case 6:
    return t8code_GPU_experiments_flow_stokes_flow_sphere_shell;
  case 7:
    return t8code_GPU_experiments_flow_around_circle_with_angular_velocity;
  default:
    SC_ABORT ("Wrong argument for flow parameter.\n");
  }
  return NULL; /* prevents compiler warning */
}

static t8code_GPU_experiments_advect_problem_t *
t8code_GPU_experiments_advect_problem_init (t8_cmesh_t cmesh, int flow, t8code_GPU_experiments_example_level_set_fn phi_0, void *ls_data,
                        int level, int maxlevel, double T, double cfl, sc_MPI_Comm comm, double band_width, int dim,
                        int dummy_op, int volume_refine)
{
  t8code_GPU_experiments_advect_problem_t *problem;
  const t8_scheme *default_scheme = t8_scheme_new_default ();
  int i;

  T8_ASSERT (1 <= dim && dim <= 3);

  /* allocate problem */
  problem = T8_ALLOC (t8code_GPU_experiments_advect_problem_t, 1);
  /* Fill problem parameters */
  problem->u = t8code_GPU_experiments_advect_choose_flow (flow); /* flow field */
  problem->flow = flow;                   /* number of the flow field */
  problem->phi_0 = phi_0;                 /* initial condition */
  problem->udata_for_phi = ls_data;       /* user data pointer passed to phi */
  problem->level = level;                 /* minimum refinement level */
  problem->maxlevel = maxlevel;           /* maximum allowed refinement level */
  problem->t = 0;                         /* start time */
  problem->T = T;                         /* end time */
  problem->delta_t = -1;                  /* delta_t, invalid value */
  problem->min_grad = 2;                  /* Coarsen an element if the gradient is smaller */
  problem->max_grad = 4;                  /* Refine an element if the gradient is larger */
  problem->min_vol = -1;                  /* Invalid start entry */
  problem->volume_refine = volume_refine; /* If greater or equal zero, refine elem only if
                                                   volume is larger than min_vol. */
  problem->cfl = cfl;                     /* cfl number  */
  problem->num_time_steps = 0;            /* current time step */
  problem->comm = comm;                   /* MPI communicator */
  problem->vtk_count = 0;                 /* number of pvtu files written */
  problem->band_width = band_width;       /* width of the refinement band around 0 level-set */
  problem->dim = dim;                     /* dimension of the mesh */
  problem->dummy_op = dummy_op;           /* If true, emulate more computational load per element */

  for (i = 0; i < ADVECT_NUM_STATS; i++) {
    sc_stats_init (&problem->stats[i], advect_stat_names[i]);
  }

  /* Construct uniform forest with ghosts */

  problem->forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  /* Initialize the element data with num_local_elements entries. */

  problem->element_data
    = t8code_GPU_experiments_advect_element_data_new (t8_forest_get_local_num_leaf_elements (problem->forest));
  problem->element_data_adapt = NULL;

  /* initialize the phi array */
  problem->phi_values
    = sc_array_new_count (sizeof (double), t8_forest_get_local_num_leaf_elements (problem->forest)
                                             + t8_forest_get_num_ghosts (problem->forest));
  problem->phi_values_adapt = NULL;
  return problem;
}

/* Project the solution at the last time step to the forest.
 * Also set the fluxes to invalid */
static void
t8code_GPU_experiments_advect_project_element_data (t8code_GPU_experiments_advect_problem_t *problem)
{
  t8_locidx_t num_local_elements, ielem;
  t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;
  int iface;

  num_local_elements = t8_forest_get_local_num_leaf_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    /* Currently the mesh does not change, thus the projected value is
     * just the computed value */
    t8code_GPU_experiments_advect_element_set_phi (problem, ielem, data->phi_new[ielem]);
    /* Set all fluxes to invalid */
    for (iface = 0; iface < data->num_faces[ielem]; iface++) {
      if (data->flux_valid[ielem][iface] >= 0) {
        /* If they are allocated (>= 0), set to allocated and not computed (=0) */
        data->flux_valid[ielem][iface] = 0;
      }
    }
  }
}

static void
t8code_GPU_experiments_advect_problem_init_elements (t8code_GPU_experiments_advect_problem_t *problem)
{
  t8_locidx_t itree, ielement, idata;
  t8_locidx_t num_trees, num_elems_in_tree;
  int iface;
  t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;
  const t8_scheme *scheme = t8_forest_get_scheme (problem->forest);
  double speed, max_speed = 0, min_diam = -1, delta_t, min_delta_t;
  t8_3D_vec u;
  double diam;
  double min_vol = 1e9;

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  /* maximum possible delta_t value */
  min_delta_t = problem->T - problem->t;
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (problem->forest, itree);
    num_elems_in_tree = t8_forest_get_tree_num_leaf_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (problem->forest, itree, ielement);
      /* Initialize the element's midpoint and volume */
      t8code_GPU_experiments_advect_compute_element_data (problem, data, idata, element, itree);
      /* Compute the minimum diameter */
      diam = t8_forest_element_diam (problem->forest, itree, element);
      T8_ASSERT (diam > 0);
      min_diam = min_diam < 0 ? diam : SC_MIN (min_diam, diam);
      /* Compute the maximum velocity */
      problem->u (data->midpoint[idata], problem->t, u);
      speed = t8_norm (u);
      max_speed = SC_MAX (max_speed, speed);

      /* Compute minimum necessary time step */
      delta_t = problem->T - problem->t;
      if (speed > 0) {
        delta_t = problem->cfl * diam / speed;
      }
      min_delta_t = SC_MIN (delta_t, min_delta_t);
      if (problem->volume_refine >= 0 && problem->min_vol <= 0) {
        /* Compute the minimum volume */
        min_vol = SC_MIN (min_vol, data->vol[idata]);
      }
      /* Set the initial condition */
      t8code_GPU_experiments_advect_element_set_phi (problem, idata, problem->phi_0 (data->midpoint[idata], 0, problem->udata_for_phi));
      /* Set the level */
      data->level[idata] = scheme->element_get_level (tree_class, element);
      /* Set the faces */
      data->num_faces[idata] = scheme->element_get_num_faces (tree_class, element);
      for (iface = 0; iface < data->num_faces[idata]; iface++) {
        /* Compute the indices of the face neighbors */
        t8code_GPU_experiments_advect_compute_face_neighbors (problem, data, idata, element, itree, iface);
        if (problem->dim > 1) {
          t8code_GPU_experiments_advect_compute_face_geometry (problem, data, idata, element, itree, iface);
        }
        data->flux_valid[idata][iface] = 0;
      }
    }
  }
  /* Exchange ghost values */
  t8_forest_ghost_exchange_data (problem->forest, problem->phi_values);

  /* Compute the timestep, this has to be done globally */
  sc_MPI_Allreduce (&min_delta_t, &problem->delta_t, 1, sc_MPI_DOUBLE, sc_MPI_MIN, problem->comm);
  if (problem->volume_refine >= 0 && problem->min_vol <= 0) {
    /* Compute the minimum volume.
     * Only in first run and only if volume refinement is active */
    sc_MPI_Allreduce (&min_vol, &problem->min_vol, 1, sc_MPI_DOUBLE, sc_MPI_MIN, problem->comm);
  }
  t8_global_essentialf ("[advect] min diam %g max flow %g  delta_t = %g\n", min_diam, max_speed, problem->delta_t);
}

static void
t8code_GPU_experiments_advect_write_vtk (t8code_GPU_experiments_advect_problem_t *problem)
{
  double *u_and_phi_array[4];
  t8_3D_vec u_temp;
  t8_locidx_t num_local_elements, ielem;
  t8_vtk_data_field_t vtk_data[5];
  const t8code_GPU_experiments_advect_element_data_t *data = problem->element_data;
  char fileprefix[BUFSIZ];
  int idim;
  double phi;

  /* Allocate num_local_elements doubles to store u and phi values */
  num_local_elements = t8_forest_get_local_num_leaf_elements (problem->forest);
  /* phi */
  u_and_phi_array[0] = T8_ALLOC_ZERO (double, num_local_elements);
  /* phi_0 */
  u_and_phi_array[1] = T8_ALLOC_ZERO (double, num_local_elements);
  /* phi - phi_0 */
  u_and_phi_array[2] = T8_ALLOC_ZERO (double, num_local_elements);
  /* u */
  u_and_phi_array[3] = T8_ALLOC_ZERO (double, 3 * (num_local_elements));

  /* Fill u and phi arrays with their values */
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    phi = t8code_GPU_experiments_advect_element_get_phi (problem, ielem);
    u_and_phi_array[0][ielem] = phi;
    u_and_phi_array[1][ielem] = problem->phi_0 (data->midpoint[ielem], problem->t, problem->udata_for_phi);
    u_and_phi_array[2][ielem] = phi - u_and_phi_array[1][ielem];
    problem->u (data->midpoint[ielem], problem->t, u_temp);
    for (idim = 0; idim < 3; idim++) {
      u_and_phi_array[3][3 * ielem + idim] = u_temp[idim];
    }
  }

  /* Write meta data for vtk */
  snprintf (vtk_data[0].description, T8_VTK_MAX_STRING_LENGTH, "Num. Solution");
  vtk_data[0].type = T8_VTK_SCALAR;
  vtk_data[0].data = u_and_phi_array[0];
  snprintf (vtk_data[1].description, T8_VTK_MAX_STRING_LENGTH, "Ana. Solution");
  vtk_data[1].type = T8_VTK_SCALAR;
  vtk_data[1].data = u_and_phi_array[1];
  snprintf (vtk_data[2].description, T8_VTK_MAX_STRING_LENGTH, "Error");
  vtk_data[2].type = T8_VTK_SCALAR;
  vtk_data[2].data = u_and_phi_array[2];
  snprintf (vtk_data[3].description, T8_VTK_MAX_STRING_LENGTH, "Flow");
  vtk_data[3].type = T8_VTK_VECTOR;
  vtk_data[3].data = u_and_phi_array[3];
  /* Write filename */
  snprintf (fileprefix, BUFSIZ, "advection_%03i", problem->vtk_count);
  /* Write vtk files */
  if (t8_forest_write_vtk_ext (problem->forest, fileprefix, 1, 1, 1, 1, 0, 0, 0, 4, vtk_data)) {
    t8_debugf ("[Advect] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Advect] Error writing to files %s\n", fileprefix);
  }
  /* clean-up */
  T8_FREE (u_and_phi_array[0]);
  T8_FREE (u_and_phi_array[1]);
  T8_FREE (u_and_phi_array[2]);
  T8_FREE (u_and_phi_array[3]);
  problem->vtk_count++;
}

#if T8_ENABLE_DEBUG
static void
t8code_GPU_experiments_advect_print_phi (t8code_GPU_experiments_advect_problem_t *problem)
{
  t8_locidx_t ielement;
  t8_locidx_t num_local_els;
  char buffer[BUFSIZ] = "";
  double phi;

  num_local_els = t8_forest_get_local_num_leaf_elements (problem->forest);
  for (ielement = 0; ielement < (t8_locidx_t) problem->phi_values->elem_count; ielement++) {
    phi = t8code_GPU_experiments_advect_element_get_phi (problem, ielement);
    snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "%.2f |%s ", phi,
              ielement == num_local_els - 1 ? "|" : "");
  }
  t8_debugf ("\t%s\n", buffer);
  /* reset buffer */
  buffer[0] = '\0';
}
#endif

static void
t8code_GPU_experiments_advect_problem_destroy (t8code_GPU_experiments_advect_problem_t **pproblem)
{
  t8code_GPU_experiments_advect_problem_t *problem;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }
  /* invalidate elements */
  t8code_GPU_experiments_advect_problem_elements_invalidate (problem);
  /* Free the element array */
  t8code_GPU_experiments_advect_element_data_destroy (&problem->element_data);
  if (problem->element_data_adapt != NULL) {
    t8code_GPU_experiments_advect_element_data_destroy (&problem->element_data_adapt);
  }
  sc_array_destroy (problem->phi_values);
  if (problem->phi_values_adapt != NULL) {
    sc_array_destroy (problem->phi_values_adapt);
  }
  /* Unref the forest */
  t8_forest_unref (&problem->forest);
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

static void
t8code_GPU_experiments_advect_solve (t8_cmesh_t cmesh, int flow, t8code_GPU_experiments_example_level_set_fn phi_0, void *ls_data,
                 const int level, const int maxlevel, double T, double cfl, sc_MPI_Comm comm, int adapt_freq,
                 int no_vtk, int vtk_freq, double band_width, int dim, int dummy_op, int volume_refine)
{
  t8code_GPU_experiments_advect_problem_t *problem;
  t8_locidx_t lelement;
  t8code_GPU_experiments_advect_element_data_t *data;
  double l_infty, L_2;
  int modulus, time_steps;
  int done = 0;
  int adapted_or_partitioned = 0;
  double total_time, solve_time = 0;
  double ghost_exchange_time, ghost_waittime, neighbor_time, flux_time;
  double vtk_time = 0;
  double start_volume, end_volume;

  /* Initialize problem */
  /* start timing */
  total_time = -sc_MPI_Wtime ();
  problem = t8code_GPU_experiments_advect_problem_init (cmesh, flow, phi_0, ls_data, level, maxlevel, T, cfl, comm, band_width, dim, dummy_op,
                                    volume_refine);
  t8code_GPU_experiments_advect_problem_init_elements (problem);

  if (maxlevel > level) {
    int ilevel;

    for (ilevel = problem->level; ilevel < problem->maxlevel; ilevel++) {
      /* initial adapt */
      t8code_GPU_experiments_advect_problem_adapt (problem, 0);
      /* repartition */
      t8code_GPU_experiments_advect_problem_partition (problem, 0);
      /* Re initialize the elements */
      t8code_GPU_experiments_advect_problem_init_elements (problem);
    }
    adapted_or_partitioned = 1;
  }
  start_volume = t8code_GPU_experiments_advect_level_set_volume (problem);
  t8_global_essentialf ("[advect] Start volume %e\n", start_volume);

#if T8_ENABLE_DEBUG
  t8code_GPU_experiments_advect_print_phi (problem);
#endif

  /* Set initialization runtime */
  sc_stats_set1 (&problem->stats[ADVECT_INIT], total_time + sc_MPI_Wtime (), advect_stat_names[ADVECT_INIT]);

  time_steps = (int) (T / problem->delta_t);
  t8_global_essentialf ("[advect] Starting with Computation. Level %i."
                        " Adaptive levels %i."
                        " End time %g. delta_t %g. cfl %g. %i time steps.\n",
                        level, maxlevel - level, T, problem->delta_t, problem->cfl, time_steps);
  T8_ASSERT (problem->delta_t > 0);
  T8_ASSERT (time_steps > 0);
  /* Controls how often we print the time step to stdout */
  modulus = SC_MAX (1, time_steps / 10);
  for (problem->num_time_steps = 0; !done; problem->num_time_steps++, problem->t += problem->delta_t) {
    if (problem->num_time_steps % modulus == modulus - 1) {
      t8_global_essentialf ("[advect] Step %i  %" T8_GLOIDX_FORMAT " elems\n", problem->num_time_steps + 1,
                            t8_forest_get_global_num_leaf_elements (problem->forest));
    }
    /* Time loop */

    /* Print vtk */
    if (!no_vtk && problem->num_time_steps % vtk_freq == 0) {
      vtk_time -= sc_MPI_Wtime ();
      t8code_GPU_experiments_advect_write_vtk (problem);
      vtk_time += sc_MPI_Wtime ();
    }
    /* Measure element count */
    sc_stats_accumulate (&problem->stats[ADVECT_ELEM_AVG], t8_forest_get_global_num_leaf_elements (problem->forest));

    solve_time -= sc_MPI_Wtime ();
    if (adapted_or_partitioned) {
      /* We changed the mesh, so that we have to compute the neighbors and the face geometry again. */
      neighbor_time = -sc_MPI_Wtime ();
      t8code_GPU_experiments_advect_compute_mesh_data (problem);
      neighbor_time += sc_MPI_Wtime ();
      sc_stats_accumulate (&problem->stats[ADVECT_NEIGHS], neighbor_time);
      /* We want to count all runs over the solver time as one */
      problem->stats[ADVECT_NEIGHS].count = 1;
    }
    /* Get the element data. This has to be done in each time step, since adapt and partition replace it.
     * Inside the element loop we only use these local variables and no problem-> or t8code calls. */
    data = problem->element_data;
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (problem->forest);
    double *phi = (double *) problem->phi_values->array;
    const int flow = problem->flow;
    const double t = problem->t;
    const double delta_t = problem->delta_t;

    /* The kernel gets a copy of the element data struct (array pointers only) as an argument.
     * Thus, it does not read the pointers from the struct in host memory, which is freed and
     * allocated again (often at the same address) in each adapt and partition. */
    t8code_GPU_experiments_advect_element_data_t kernel_data = *data;

    flux_time = -sc_MPI_Wtime ();
#pragma omp target teams distribute parallel for firstprivate (kernel_data)
    for (lelement = 0; lelement < num_local_elements; lelement++) {
      /* element loop */
      /* Compute the fluxes across all faces */
      t8code_GPU_experiments_advect_compute_fluxes (&kernel_data, lelement, phi, dim, flow, t);
      if (dummy_op) {
        /* simulate more load per element */
        t8code_GPU_experiments_advect_dummy_op (&kernel_data, lelement);
      }
      /* Compute time step for this element */
      t8code_GPU_experiments_advect_advance_element (&kernel_data, lelement, phi, delta_t);
    } // end of element loop
    flux_time += sc_MPI_Wtime ();
    sc_stats_accumulate (&problem->stats[ADVECT_FLUX], flux_time);
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_FLUX].count = 1;
    adapted_or_partitioned = 0;
    /* Store the advanced phi value in each element */
    t8code_GPU_experiments_advect_project_element_data (problem);
    solve_time += sc_MPI_Wtime ();
    if (maxlevel > level) {
      /* Adapt the mesh after adapt_freq time steps */
      if (problem->num_time_steps % adapt_freq == adapt_freq - 1) {
        adapted_or_partitioned = 1;
        t8code_GPU_experiments_advect_problem_adapt (problem, 1);
        t8code_GPU_experiments_advect_problem_partition (problem, 1);
      }
    }

    /* Exchange ghost values */
    ghost_exchange_time = -sc_MPI_Wtime ();
    t8_forest_ghost_exchange_data (problem->forest, problem->phi_values);
    ghost_exchange_time += sc_MPI_Wtime ();
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_EXCHANGE], ghost_exchange_time);
    ghost_waittime = t8_forest_profile_get_ghostexchange_waittime (problem->forest);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_WAIT], ghost_waittime);
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_GHOST_EXCHANGE].count = 1;
    problem->stats[ADVECT_GHOST_WAIT].count = 1;

    if (problem->t + problem->delta_t > problem->T) {
      /* Ensure that the last time step is always the given end time */
      problem->delta_t = problem->T - problem->t;
    }
    /* Check whether we are finished */
    if (problem->t >= problem->T) {
      done = 1;
    }
  } /* End element loop */
  if (!no_vtk) {
    vtk_time -= sc_MPI_Wtime ();
    /* Print last time step vtk */
    t8code_GPU_experiments_advect_write_vtk (problem);
    vtk_time += sc_MPI_Wtime ();
  }
  /* Compute runtime */
  total_time += sc_MPI_Wtime ();
  sc_stats_set1 (&problem->stats[ADVECT_TOTAL], total_time, advect_stat_names[ADVECT_TOTAL]);
  sc_stats_set1 (&problem->stats[ADVECT_SOLVE], solve_time, advect_stat_names[ADVECT_SOLVE]);
  sc_stats_set1 (&problem->stats[ADVECT_IO], vtk_time, advect_stat_names[ADVECT_IO]);
  /* Compute volume loss */

  end_volume = t8code_GPU_experiments_advect_level_set_volume (problem);
  t8_global_essentialf ("[advect] End volume %e\n", end_volume);

  sc_stats_set1 (&problem->stats[ADVECT_VOL_LOSS], 100 * (1 - end_volume / start_volume),
                 advect_stat_names[ADVECT_VOL_LOSS]);

  /* Compute l_infty error */
  l_infty = t8code_GPU_experiments_advect_l_infty_rel (problem, phi_0, 0.025);
  L_2 = t8code_GPU_experiments_advect_l_2_rel (problem, phi_0, 0.025);
  t8_global_essentialf ("[advect] Done. t = %g \t l_infty error:\t%e\tL_2:\t%e\n", problem->t, l_infty, L_2);

  sc_stats_set1 (&problem->stats[ADVECT_ERROR_INF], l_infty, advect_stat_names[ADVECT_ERROR_INF]);
  sc_stats_set1 (&problem->stats[ADVECT_ERROR_2], L_2, advect_stat_names[ADVECT_ERROR_2]);
  sc_stats_compute (problem->comm, ADVECT_NUM_STATS, problem->stats);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, ADVECT_NUM_STATS, problem->stats, 1, 1);
  /* clean-up */
  t8code_GPU_experiments_advect_problem_destroy (&problem);
}

int
main (int argc, char *argv[])
{
  int mpiret;
  sc_options_t *opt;
  char help[BUFSIZ];
  const char *mshfile = NULL;
  int level, reflevel, dim, cube_type, dummy_op;
  int parsed, helpme, no_vtk, vtk_freq, adapt_freq;
  int volume_refine;
  int flow_arg, use_cad_geometry;
  double T, cfl, band_width;
  t8code_GPU_experiments_levelset_sphere_data_t ls_data;
  /* brief help message */

  /* long help message */

  snprintf (help, BUFSIZ,
            "This program solves the advection equation on "
            "a given geometry.\n");
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
#if T8_ENABLE_DEBUG
  t8_init (SC_LP_DEBUG);
#else
  t8_init (SC_LP_ESSENTIAL);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'u', "flow", &flow_arg, 0,
                      "Choose the flow field u.\n"
                      "\t\t1 - Constant 1 in x-direction.\n"
                      "\t\t2 - Constant 1 in x,y, and z.\n"
                      "\t\t3 - A turbulent flow in a cube with zero outflow.\n"
                      "\t\t\tIt reverses direction at t = 0.5.\n"
                      "\t\t4 - 2D rotation around (0.5,0.5).\n"
                      "\t\t5 - 2D flow around circle at (0.5,0.5)"
                      "with radius 0.15.\n"
                      "\t\t6 - A solution to the stokes equation on a spherical shell.\n"
                      "\t\t7 - Flow past a rotating cylinder of radius of 0.5"
                      " around the z-axis. This flow is defined for a specific"
                      " mesh, which can be generated with Gmsh and the .geo"
                      " files 't8_advection_generate_channel.geo' and"
                      " 't8_advection_generate_channel_2d.geo'. These meshes"
                      " can also be used with the curved geometry.\n");
  sc_options_add_int (opt, 'l', "level", &level, 0, "The minimum refinement level of the mesh.");
  sc_options_add_int (opt, 'r', "rlevel", &reflevel, 0, "The number of adaptive refinement levels.");
  sc_options_add_int (opt, 'e', "elements", &cube_type, -1,
                      "If specified the coarse mesh is a hypercube\n\t\t\t\t     consisting of the"
                      " following elements:\n"
                      "\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n"
                      "\t\t7 - triangle/quad (hybrid 2d).\n"
                      "\t\t8 - tet/hex/prism (hybrid 3d).\n"
                      "\t\t9 - pyramid.\n");
  sc_options_add_string (opt, 'f', "mshfile", &mshfile, NULL,
                         "If specified, the cmesh is constructed from a .msh file with "
                         "the given prefix.\n\t\t\t\t     The files must end in .msh "
                         "and be in ASCII format version 2. -d must be specified.");
  sc_options_add_int (opt, 'd', "dim", &dim, -1, "In combination with -f: The dimension of the mesh. 1 <= d <= 3.");

  sc_options_add_switch (opt, 'c', "cad", &use_cad_geometry,
                         "In combination with -f: Use the cad geometry, only viable if a "
                         ".brep file of the same name is present.");

  sc_options_add_double (opt, 'T', "end-time", &T, 1, "The duration of the simulation. Default: 1");

  sc_options_add_double (opt, 'C', "CFL", &cfl, 1, "The cfl number to use. Default: 1");
  sc_options_add_double (opt, 'b', "band-width", &band_width, 1,
                         "Control the width of the refinement band around\n"
                         "\t\t\t\t     the zero level-set. Default 1.");

  sc_options_add_int (opt, 'a', "adapt-freq", &adapt_freq, 1,
                      "Controls how often the mesh is readapted. "
                      "A value of i means, every i-th time step.");

  sc_options_add_int (opt, 'v', "vtk-freq", &vtk_freq, 1,
                      "How often the vtk output is produced "
                      "\n\t\t\t\t     (after how many time steps). "
                      "A value of 0 is equivalent to using -o.");

  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk,
                         "Suppress vtk output. "
                         "Overwrites any -v setting.");

  sc_options_add_switch (opt, 's', "simulate", &dummy_op,
                         "Simulate more load per element. "
                         "In each iteration, useless dummy operations\n "
                         "\t\t\t\t     are performed per element. Decreases the "
                         "performance!");
  sc_options_add_double (opt, 'X', "Xcoord", &ls_data.M[0], 0.6,
                         "The X-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'Y', "Ycoord", &ls_data.M[1], 0.6,
                         "The Y-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'Z', "Zcoord", &ls_data.M[2], 0.6,
                         "The Z-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'R', "Radius", &ls_data.radius, 0.25,
                         "The radius of the Sphere."
                         "Default is 0.25.");

  sc_options_add_int (opt, 'V', "volume-refine", &volume_refine, -1,
                      "Refine elements close to the 0 level-set only "
                      "if their volume is smaller than the l+V-times refined\n"
                      "\t\t\t\t     smallest element int the mesh.");

  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 1 <= flow_arg && flow_arg <= 7 && 0 <= level && 0 <= reflevel && 0 <= vtk_freq
           && ((mshfile != NULL && 0 < dim && dim <= 3) || (1 <= cube_type && cube_type <= 9)) && band_width >= 0) {
    t8_cmesh_t cmesh;
    if (mshfile == NULL) {
      switch (cube_type) {
      case 7:
        dim = 2;
        break;
      case 8:
        dim = 3;
        break;
      case 9:
        dim = 3;
        break;
      default:
        dim = t8_eclass_to_dimension[cube_type];
        T8_ASSERT (cube_type < 7);
      }
    }
    /* Set level-set midpoint coordinates to zero for unused dimensions. */
    if (cube_type == 2 || cube_type == 3 || cube_type == 7) {
      ls_data.M[2] = 0;
    }
    if (cube_type == 1) {
      ls_data.M[1] = ls_data.M[2] = 0;
    }

    cmesh = t8code_GPU_experiments_advect_create_cmesh (sc_MPI_COMM_WORLD, cube_type, mshfile, level, dim, use_cad_geometry);
    if (!no_vtk) {
      t8_cmesh_vtk_write_file (cmesh, "advection_cmesh");
    }
    /* Computation */
    t8code_GPU_experiments_advect_solve (cmesh, flow_arg, t8code_GPU_experiments_levelset_sphere, &ls_data, level, level + reflevel, T, cfl, sc_MPI_COMM_WORLD,
                     adapt_freq, no_vtk, vtk_freq, band_width, dim, dummy_op, volume_refine);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR:Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
