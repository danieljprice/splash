/*
 * HDF5 read utilities for Phantom / sphNG particle data (C API only).
 *
 * Standalone Phantom HDF5 reader and chemistry sidecar support for binary
 * Phantom dumps. Called from read_data_phantom_hdf5.f90 and
 * read_composition_hdf5.f90 via iso_c_binding.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <hdf5.h>
#include "hdf5_helper_utils.h"

#define PHANTOM_NAMELEN 256
#define PHANTOM_PARTICLES "particles"
#define PHANTOM_HEADER "header"

/* Fortran callbacks */
void set_blocklabel_phantom(int *icol, int *irank, char *name);
void read_phantom_hdf5_data_fromc(int *icol, int *npart, double temparr[], char *name);
void set_extra_column_label_phantom(int *icol, char *name);
void read_extra_column_fromc_phantom(int *icol, int *npart, double temparr[], int *icomp_col_start);

/* bind(c) entry points */
int phantom_hdf5_is_phantom_file(char *filename);
void phantom_hdf5_header(char *filename, double *time, int *npart, int *ncol, int *has_header,
                         int *ndim, int *ierr);
void phantom_hdf5_data(char *filename, int *npart, int *ncol, int *isrequired, int *ierr);
void phantom_hdf5_extra_check(char *filename, int *ntotal, int *ncomp, int *npart_file, int *ierr);
void phantom_hdf5_read_particle_ids(char *filename, int *npart, int *ids, int *ierr);
void phantom_hdf5_extra_read(char *filename, int *ntotal, int *icomp_col_start, int *ncomp,
                             int *isrequired, int *ierr);
static const herr_t HDF5_error = -1;
static int extra_nonfinite_nvals = 0;
static unsigned char *extra_particle_bad = NULL;
static int extra_particle_bad_np = 0;
static int extra_icomp_col_start = 0;

/*
 * allocate or reset per-particle NaN/Inf tracking for sidecar reads
 */
static void reset_extra_nonfinite_tracking(int np)
{
    if (np != extra_particle_bad_np)
    {
        free(extra_particle_bad);
        extra_particle_bad = NULL;
        extra_particle_bad_np = 0;
    }
    if (np > 0 && !extra_particle_bad)
    {
        extra_particle_bad = (unsigned char *)calloc((size_t)np, sizeof(unsigned char));
        if (extra_particle_bad)
            extra_particle_bad_np = np;
    }
    extra_nonfinite_nvals = 0;
}

/*
 * replace NaN/Inf in a sidecar column buffer with zero and track bad particles
 */
static void sanitize_extra_buffer(double *buf, int n)
{
    int i;

    if (n > 0 && n != extra_particle_bad_np)
        reset_extra_nonfinite_tracking(n);

    for (i = 0; i < n; i++)
    {
        if (!isfinite(buf[i]))
        {
            buf[i] = 0.0;
            extra_nonfinite_nvals++;
            if (extra_particle_bad && i < extra_particle_bad_np)
                extra_particle_bad[i] = 1;
        }
    }
}

/*
 * print summary warning if sidecar chemistry columns contained NaN/Inf values
 */
static void print_extra_nonfinite_warning(void)
{
    int i, nbad = 0;

    if (!extra_particle_bad || extra_nonfinite_nvals <= 0)
        return;
    for (i = 0; i < extra_particle_bad_np; i++)
        if (extra_particle_bad[i])
            nbad++;
    if (nbad > 0)
        printf(" WARNING: chemistry sidecar had NaN/Inf on %i particle%s (%i values across"
               " species; set to 0)\n",
               nbad, nbad == 1 ? "" : "s", extra_nonfinite_nvals);
    free(extra_particle_bad);
    extra_particle_bad = NULL;
    extra_particle_bad_np = 0;
    extra_nonfinite_nvals = 0;
}

/*
 * open an HDF5 group (HDF5 1.8+ compatible wrapper)
 */
static hid_t open_group(hid_t loc_id, const char *name)
{
#if H5_VERSION_GE(1, 8, 0)
    return H5Gopen2(loc_id, name, H5P_DEFAULT);
#else
    return H5Gopen(loc_id, name);
#endif
}

/*
 * open an HDF5 dataset (HDF5 1.8+ compatible wrapper)
 */
static hid_t open_dataset(hid_t loc_id, const char *name)
{
#if H5_VERSION_GE(1, 8, 0)
    return H5Dopen2(loc_id, name, H5P_DEFAULT);
#else
    return H5Dopen(loc_id, name);
#endif
}

/*
 * convert string to lower case in place
 */
static void strtolower(char *s)
{
    int i;
    for (i = 0; s[i] != '\0'; i++)
        s[i] = (char)tolower((unsigned char)s[i]);
}

/*
 * trim leading and trailing whitespace from a string in place
 */
static void trim_name(char *s)
{
    int start = 0, end;
    if (!s)
        return;
    while (s[start] != '\0' && isspace((unsigned char)s[start]))
        start++;
    if (start > 0)
        memmove(s, s + start, strlen(s + start) + 1);
    end = (int)strlen(s) - 1;
    while (end >= 0 && isspace((unsigned char)s[end]))
        s[end--] = '\0';
}

/*
 * return true if dataset name is a core hydro field (not a chemistry sidecar column)
 */
static int name_in_skip_list(const char *name)
{
    static const char *skip[] = {
        "x", "y", "z", "vx", "vy", "vz", "vxyz", "vels", "rho", "pmass", "mass",
        "time", "id", "iorig", "iphase", "tag", "av", "u", "divv", "dt", "alpha", "beta",
        "temp", "temperature",
        "placeholder", /* phantom dummy column (often all NaN) */
        NULL
    };
    char lname[PHANTOM_NAMELEN];
    int i;

    strncpy(lname, name, PHANTOM_NAMELEN - 1);
    lname[PHANTOM_NAMELEN - 1] = '\0';
    trim_name(lname);
    strtolower(lname);
    for (i = 0; skip[i] != NULL; i++)
        if (strcmp(lname, skip[i]) == 0)
            return 1;
    return 0;
}

/*
 * return number of particles represented by a dataset dataspace (1D or 3xN / Nx3)
 */
static int dataset_length(hid_t dataspace_id, int *rank_out, hsize_t dims_out[2])
{
    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    hsize_t dims[2] = {0, 0}, maxdims[2] = {0, 0};

    if (rank < 1)
        return 0;
    H5Sget_simple_extent_dims(dataspace_id, dims, maxdims);
    dims_out[0] = dims[0];
    dims_out[1] = (rank > 1) ? dims[1] : 0;
    *rank_out = rank;
    if (rank == 1)
        return (int)dims[0];
    if (rank == 2)
    {
        if (dims[0] == 3)
            return (int)dims[1];
        if (dims[1] == 3)
            return (int)dims[0];
    }
    return 0;
}

/*
 * return number of splash columns occupied by a dataset of given rank/shape
 */
static int ncol_for_rank(int rank, hsize_t d0, hsize_t d1)
{
    if (rank == 1)
        return 1;
    if (rank == 2 && (d0 == 3 || d1 == 3))
        return 3;
    return 0;
}

/*
 * read a 1D particle dataset into a double buffer (handles float/int types)
 */
static int read_scalar_dataset(hid_t dataset_id, hid_t type_id, int npart, double *buffer)
{
    int i;
    herr_t status;

    if (H5Tequal(type_id, H5T_NATIVE_DOUBLE))
    {
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
        return (status == HDF5_error);
    }
    if (H5Tequal(type_id, H5T_NATIVE_FLOAT))
    {
        float *tmp = (float *)malloc((size_t)npart * sizeof(float));
        if (!tmp)
            return 1;
        status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
        if (status != HDF5_error)
            for (i = 0; i < npart; i++)
                buffer[i] = (double)tmp[i];
        free(tmp);
        return (status == HDF5_error);
    }
    if (H5Tequal(type_id, H5T_NATIVE_INT) || H5Tequal(type_id, H5T_STD_I32LE))
    {
        int *tmp = (int *)malloc((size_t)npart * sizeof(int));
        if (!tmp)
            return 1;
        status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
        if (status != HDF5_error)
            for (i = 0; i < npart; i++)
                buffer[i] = (double)tmp[i];
        free(tmp);
        return (status == HDF5_error);
    }
    if (H5Tequal(type_id, H5T_NATIVE_LLONG) || H5Tequal(type_id, H5T_STD_I64LE))
    {
        long long *tmp = (long long *)malloc((size_t)npart * sizeof(long long));
        if (!tmp)
            return 1;
        status = H5Dread(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
        if (status != HDF5_error)
            for (i = 0; i < npart; i++)
                buffer[i] = (double)tmp[i];
        free(tmp);
        return (status == HDF5_error);
    }
    return 1;
}

/*
 * read one component of a 3xN or Nx3 vector dataset into a double buffer
 */
static int read_vector_component(hid_t dataset_id, hid_t type_id, int npart, int icomp,
                                 hsize_t d0, double *buffer)
{
    hid_t filespace, memspace;
    hsize_t offset[2], count[2];
    herr_t status;
    int i;

    filespace = H5Dget_space(dataset_id);
    if (d0 == 3)
    {
        offset[0] = (hsize_t)icomp;
        offset[1] = 0;
        count[0] = 1;
        count[1] = (hsize_t)npart;
    }
    else
    {
        offset[0] = 0;
        offset[1] = (hsize_t)icomp;
        count[0] = (hsize_t)npart;
        count[1] = 1;
    }
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    {
        hsize_t mem_n = (d0 == 3) ? count[1] : count[0];
        memspace = H5Screate_simple(1, &mem_n, NULL);
    }

    if (H5Tequal(type_id, H5T_NATIVE_DOUBLE))
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, buffer);
    else if (H5Tequal(type_id, H5T_NATIVE_FLOAT))
    {
        float *tmp = (float *)malloc((size_t)npart * sizeof(float));
        status = HDF5_error;
        if (tmp)
        {
            status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmp);
            if (status != HDF5_error)
                for (i = 0; i < npart; i++)
                    buffer[i] = (double)tmp[i];
            free(tmp);
        }
    }
    else
        status = HDF5_error;

    H5Sclose(memspace);
    H5Sclose(filespace);
    return (status == HDF5_error);
}

/*
 * get number of particles from the x dataset in a particles group
 */
static int npart_from_x(hid_t particles_id)
{
    hid_t dataset_id, dataspace_id;
    int rank, npart;
    hsize_t dims[2];

    if (!checkfordataset(particles_id, "x"))
        return 0;
    dataset_id = open_dataset(particles_id, "x");
    if (dataset_id == HDF5_error)
        return 0;
    dataspace_id = H5Dget_space(dataset_id);
    npart = dataset_length(dataspace_id, &rank, dims);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return npart;
}

/* metadata for one dataset in the particles group (used for sorting/iteration) */
typedef struct
{
    char name[PHANTOM_NAMELEN];
    int rank;
    hsize_t dim0, dim1;
    int ncols;
    int file_index;
} phantom_dataset_t;

/*
 * copy dataset name to lower case in a fixed-size buffer
 */
static void copy_lower(const char *name, char *lname)
{
    int i;

    strncpy(lname, name, PHANTOM_NAMELEN - 1);
    lname[PHANTOM_NAMELEN - 1] = '\0';
    for (i = 0; lname[i] != '\0'; i++)
        lname[i] = (char)tolower((unsigned char)lname[i]);
}

/*
 * sort key for standalone phantom_hdf5 read: xyz, hydro, av, time/thermal, abundances
 */
static int standalone_column_sort_key(const char *name)
{
    static const char *xyz[] = {"x", "y", "z", NULL};
    static const char *hydro[] = {
        "vx", "vy", "vz", "vxyz", "vels", "rho", "density", "pmass", "mass",
        "u", "divv", "dt", "alpha", "beta", "id", "iorig", "iphase", "tag", NULL
    };
    static const char *time_thermal[] = {"time", "temperature", "temp", NULL};
    char lname[PHANTOM_NAMELEN];
    int i;

    copy_lower(name, lname);
    for (i = 0; xyz[i] != NULL; i++)
        if (strcmp(lname, xyz[i]) == 0)
            return i;
    for (i = 0; hydro[i] != NULL; i++)
        if (strcmp(lname, hydro[i]) == 0)
            return 100 + i;
    if (strcmp(lname, "av") == 0)
        return 1000;
    for (i = 0; time_thermal[i] != NULL; i++)
        if (strcmp(lname, time_thermal[i]) == 0)
            return 1100 + i;
    return 2000;
}

static int sort_standalone_columns;

/*
 * qsort comparator for phantom_dataset_t (file order or standalone sort key)
 */
static int compare_phantom_datasets(const void *a, const void *b)
{
    const phantom_dataset_t *da = (const phantom_dataset_t *)a;
    const phantom_dataset_t *db = (const phantom_dataset_t *)b;
    int ka, kb;

    if (!sort_standalone_columns)
        return da->file_index - db->file_index;

    ka = standalone_column_sort_key(da->name);
    kb = standalone_column_sort_key(db->name);
    if (ka != kb)
        return ka - kb;
    return strcmp(da->name, db->name);
}

/*
 * read or count one dataset in the particles group; callback to Fortran for data/labels
 */
static int process_one_particle_dataset(hid_t particles_id, const phantom_dataset_t *ds,
                                        int skip_core, int count_only, int *icol, int np,
                                        int *ncols, int isrequired[])
{
    hid_t dataset_id, dataspace_id, type_id;
    char name[PHANTOM_NAMELEN];
    hsize_t dims[2];
    int rank, nlen, k;
    double *buffer = NULL;

    strncpy(name, ds->name, PHANTOM_NAMELEN - 1);
    name[PHANTOM_NAMELEN - 1] = '\0';
    rank = ds->rank;
    dims[0] = ds->dim0;
    dims[1] = ds->dim1;

    dataset_id = open_dataset(particles_id, name);
    if (dataset_id == HDF5_error)
        return 0;
    dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0)
    {
        H5Dclose(dataset_id);
        return 0;
    }
    type_id = H5Dget_type(dataset_id);
    if (type_id < 0)
    {
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return 0;
    }
    nlen = dataset_length(dataspace_id, &rank, dims);
    if (nlen <= 0 || (np > 0 && nlen != np))
    {
        H5Tclose(type_id);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return 0;
    }

    if (rank == 1)
    {
        (*icol)++;
        if (skip_core)
            set_extra_column_label_phantom(icol, name);
        else if (!count_only)
        {
            int irank = 1;
            set_blocklabel_phantom(icol, &irank, name);
        }
        if (!count_only)
        {
            int req = (isrequired == NULL || isrequired[*icol - 1]);
            if (req)
            {
                buffer = (double *)malloc((size_t)np * sizeof(double));
                if (!buffer)
                {
                    H5Tclose(type_id);
                    H5Sclose(dataspace_id);
                    H5Dclose(dataset_id);
                    return 1;
                }
                if (read_scalar_dataset(dataset_id, type_id, np, buffer) == 0)
                {
                    if (skip_core)
                        sanitize_extra_buffer(buffer, np);
                    if (skip_core)
                        read_extra_column_fromc_phantom(icol, &np, buffer, &extra_icomp_col_start);
                    else
                        read_phantom_hdf5_data_fromc(icol, &np, buffer, name);
                }
                free(buffer);
                buffer = NULL;
            }
        }
    }
    else
    {
        int irank = 3;
        static const char comp_xyz[] = "xyz";
        if (!count_only)
            set_blocklabel_phantom(icol, &irank, name);
        for (k = 0; k < 3; k++)
        {
            (*icol)++;
            if (skip_core)
            {
                char lname[PHANTOM_NAMELEN];
                snprintf(lname, PHANTOM_NAMELEN, "%s_%c", name, comp_xyz[k]);
                set_extra_column_label_phantom(icol, lname);
            }
            if (!count_only)
            {
                int req = (isrequired == NULL || isrequired[*icol - 1]);
                if (req)
                {
                    if (!buffer)
                        buffer = (double *)malloc((size_t)np * sizeof(double));
                    if (buffer && read_vector_component(dataset_id, type_id, np, k, ds->dim0, buffer) == 0)
                    {
                        if (skip_core)
                            sanitize_extra_buffer(buffer, np);
                        if (skip_core)
                            read_extra_column_fromc_phantom(icol, &np, buffer, &extra_icomp_col_start);
                        else
                            read_phantom_hdf5_data_fromc(icol, &np, buffer, name);
                    }
                }
            }
        }
        if (buffer)
        {
            free(buffer);
            buffer = NULL;
        }
    }

    *ncols += ds->ncols;
    H5Tclose(type_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return 0;
}

/*
 * walk all datasets in the particles group; count or read columns
 */
static int process_particles_group(hid_t particles_id, int skip_core, int count_only,
                                   int *ncol, int *npart, int isrequired[])
{
    hsize_t nobj;
    int i, icol = 0, ncols = 0, np = 0;
    char name[PHANTOM_NAMELEN];
    phantom_dataset_t *datasets = NULL;
    int ndatasets = 0, ndatasets_alloc = 0;

    if (H5Gget_num_objs(particles_id, &nobj) == HDF5_error)
        return 1;

    if (*npart > 0)
        np = *npart;
    else
        np = npart_from_x(particles_id);

    for (i = 0; i < (int)nobj; i++)
    {
        hid_t dataset_id, dataspace_id, type_id;
        hsize_t dims[2];
        int rank, nlen, nadd;

        if (H5Gget_objname_by_idx(particles_id, (hsize_t)i, name, PHANTOM_NAMELEN) < 0)
            continue;
        name[PHANTOM_NAMELEN - 1] = '\0';
        trim_name(name);
        if (H5Gget_objtype_by_idx(particles_id, (hsize_t)i) != H5G_DATASET)
            continue;
        if (skip_core && name_in_skip_list(name))
            continue;

        dataset_id = open_dataset(particles_id, name);
        if (dataset_id == HDF5_error)
            continue;
        dataspace_id = H5Dget_space(dataset_id);
        if (dataspace_id < 0)
        {
            H5Dclose(dataset_id);
            continue;
        }
        type_id = H5Dget_type(dataset_id);
        if (type_id < 0)
        {
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            continue;
        }
        nlen = dataset_length(dataspace_id, &rank, dims);
        nadd = ncol_for_rank(rank, dims[0], dims[1]);

        if (nlen <= 0 || nadd == 0 || (np > 0 && nlen != np))
        {
            H5Tclose(type_id);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            continue;
        }
        if (np == 0)
            np = nlen;

        if (ndatasets >= ndatasets_alloc)
        {
            int nnew = (ndatasets_alloc == 0) ? 64 : ndatasets_alloc * 2;
            phantom_dataset_t *tmp =
                (phantom_dataset_t *)realloc(datasets, (size_t)nnew * sizeof(phantom_dataset_t));
            if (!tmp)
            {
                H5Tclose(type_id);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                free(datasets);
                return 1;
            }
            datasets = tmp;
            ndatasets_alloc = nnew;
        }
        strncpy(datasets[ndatasets].name, name, PHANTOM_NAMELEN - 1);
        datasets[ndatasets].name[PHANTOM_NAMELEN - 1] = '\0';
        datasets[ndatasets].rank = rank;
        datasets[ndatasets].dim0 = dims[0];
        datasets[ndatasets].dim1 = dims[1];
        datasets[ndatasets].ncols = nadd;
        datasets[ndatasets].file_index = i;
        ndatasets++;

        H5Tclose(type_id);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
    }

    sort_standalone_columns = (skip_core == 0);
    if (ndatasets > 1)
        qsort(datasets, (size_t)ndatasets, sizeof(phantom_dataset_t), compare_phantom_datasets);

    for (i = 0; i < ndatasets; i++)
    {
        if (process_one_particle_dataset(particles_id, &datasets[i], skip_core, count_only, &icol,
                                         np, &ncols, isrequired) != 0)
        {
            free(datasets);
            return 1;
        }
    }
    free(datasets);

    if (ncol)
        *ncol = icol; /* labelled column count (must match set_extra_column_label calls) */
    if (npart)
        *npart = np;
    return 0;
}

/*
 * read element 0 from a 1D (or scalar) dataset into *value
 */
static int read_first_element(hid_t loc_id, const char *dsetname, double *value)
{
    hid_t dataset_id, filespace, memspace;
    hsize_t offset[1] = {0}, count[1] = {1};
    herr_t status;
    int rank;
    double dval;
    float fval;

    *value = 0.0;
    if (!checkfordataset(loc_id, (char *)dsetname))
        return 0;

    dataset_id = open_dataset(loc_id, dsetname);
    if (dataset_id == HDF5_error)
        return 1;

    filespace = H5Dget_space(dataset_id);
    rank = H5Sget_simple_extent_ndims(filespace);
    if (rank < 1)
    {
        H5Sclose(filespace);
        H5Dclose(dataset_id);
        return 1;
    }

    memspace = H5Screate_simple(1, count, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &dval);
    if (status != HDF5_error)
        *value = dval;
    else
    {
        status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, &fval);
        if (status != HDF5_error)
            *value = (double)fval;
    }

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset_id);
    return (status == HDF5_error);
}

/*
 * read simulation time from the header/time dataset
 */
static int read_header_time(hid_t header_id, double *time)
{
    return read_first_element(header_id, "time", time);
}

/*
 * return 1 if filename looks like Phantom/sphNG HDF5 (not Gadget Header group)
 */
int phantom_hdf5_is_phantom_file(char *filename)
{
    hid_t file_id, particles_id;
    int is_phantom = 0;

    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == HDF5_error)
        return 0;
    if (checkfordataset(file_id, "Header"))
    {
        H5Fclose(file_id);
        return 0;
    }
    if (checkfordataset(file_id, PHANTOM_HEADER))
        is_phantom = 1;
    if (checkfordataset(file_id, PHANTOM_PARTICLES))
    {
        particles_id = open_group(file_id, PHANTOM_PARTICLES);
        if (particles_id != HDF5_error)
        {
            if (npart_from_x(particles_id) > 0)
                is_phantom = 1;
            H5Gclose(particles_id);
        }
    }
    H5Fclose(file_id);
    return is_phantom;
}

/*
 * read Phantom HDF5 header: time, npart, ncol (count-only pass over particles group)
 */
void phantom_hdf5_header(char *filename, double *time, int *npart, int *ncol, int *has_header,
                         int *ndim, int *ierr)
{
    hid_t file_id, header_id, particles_id;
    int np = 0, nc = 0;

    *ierr = 0;
    *time = 0.0;
    *ndim = 3;
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == HDF5_error)
    {
        *ierr = 1;
        return;
    }
    *has_header = checkfordataset(file_id, PHANTOM_HEADER);
    if (*has_header)
    {
        header_id = open_group(file_id, PHANTOM_HEADER);
        if (header_id != HDF5_error)
        {
            read_header_time(header_id, time);
            H5Gclose(header_id);
        }
    }
    if (!checkfordataset(file_id, PHANTOM_PARTICLES))
    {
        printf(" ERROR: \"%s\" group not found in Phantom HDF5 file\n", PHANTOM_PARTICLES);
        H5Fclose(file_id);
        *ierr = 2;
        return;
    }
    particles_id = open_group(file_id, PHANTOM_PARTICLES);
    if (particles_id == HDF5_error)
    {
        H5Fclose(file_id);
        *ierr = 2;
        return;
    }
    process_particles_group(particles_id, 0, 1, &nc, &np, NULL);
    if (np <= 0)
    {
        printf(" ERROR: x positions not found in Phantom HDF5 file\n");
        H5Gclose(particles_id);
        H5Fclose(file_id);
        *ierr = 3;
        return;
    }
    if (*time == 0.0)
        read_first_element(particles_id, "time", time);
    *npart = np;
    *ncol = nc;
    H5Gclose(particles_id);
    H5Fclose(file_id);
}

/*
 * read all required columns from a standalone Phantom HDF5 dump into dat
 */
void phantom_hdf5_data(char *filename, int *npart, int *ncol, int *isrequired, int *ierr)
{
    hid_t file_id, particles_id;

    *ierr = 0;
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == HDF5_error)
    {
        *ierr = 1;
        return;
    }
    particles_id = open_group(file_id, PHANTOM_PARTICLES);
    if (particles_id == HDF5_error)
    {
        H5Fclose(file_id);
        *ierr = 2;
        return;
    }
    process_particles_group(particles_id, 0, 0, ncol, npart, isrequired);
    H5Gclose(particles_id);
    H5Fclose(file_id);
}

/*
 * check chemistry sidecar: count extra species columns and npart in dump_XXXXX.h5
 */
void phantom_hdf5_extra_check(char *filename, int *ntotal, int *ncomp, int *npart_file, int *ierr)
{
    hid_t file_id, particles_id;
    int nc = 0, np = 0;

    *ierr = 0;
    *ncomp = 0;
    *npart_file = 0;
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == HDF5_error)
    {
        *ierr = 1;
        return;
    }
    if (!checkfordataset(file_id, PHANTOM_PARTICLES))
    {
        H5Fclose(file_id);
        *ierr = 2;
        return;
    }
    particles_id = open_group(file_id, PHANTOM_PARTICLES);
    if (particles_id == HDF5_error)
    {
        H5Fclose(file_id);
        *ierr = 2;
        return;
    }
    process_particles_group(particles_id, 1, 1, &nc, &np, NULL);
    *ncomp = nc;
    *npart_file = np;
    if (np > 0 && *ntotal > 0 && np > *ntotal)
    {
        printf(" ERROR: Phantom HDF5 sidecar npart (%i) > dump npart (%i)\n", np, *ntotal);
        *ierr = 3;
    }
    else if (np > 0 && *ntotal > 0 && np < *ntotal)
    {
        printf(" WARNING: Phantom HDF5 sidecar has %i particles; dump has %i slots", np, *ntotal);
        printf(" (dead/accreted slots in dump; matching sidecar by particle id)\n");
    }
    H5Gclose(particles_id);
    H5Fclose(file_id);
}

/*
 * read particle ids from sidecar for id-based scatter onto binary dump rows
 */
void phantom_hdf5_read_particle_ids(char *filename, int *npart, int *ids, int *ierr)
{
    hid_t file_id, particles_id, dataset_id, dataspace_id, type_id;
    int rank, nlen, i;
    hsize_t dims[2];

    *ierr = 0;
    if (*npart <= 0 || ids == NULL)
        return;

    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == HDF5_error)
    {
        *ierr = 1;
        return;
    }
    particles_id = open_group(file_id, PHANTOM_PARTICLES);
    if (particles_id == HDF5_error)
    {
        H5Fclose(file_id);
        *ierr = 2;
        return;
    }
    if (!checkfordataset(particles_id, "id"))
    {
        H5Gclose(particles_id);
        H5Fclose(file_id);
        *ierr = 3;
        return;
    }
    dataset_id = open_dataset(particles_id, "id");
    if (dataset_id == HDF5_error)
    {
        H5Gclose(particles_id);
        H5Fclose(file_id);
        *ierr = 3;
        return;
    }
    dataspace_id = H5Dget_space(dataset_id);
    type_id = H5Dget_type(dataset_id);
    nlen = dataset_length(dataspace_id, &rank, dims);
    if (nlen != *npart || rank != 1)
    {
        H5Tclose(type_id);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Gclose(particles_id);
        H5Fclose(file_id);
        *ierr = 4;
        return;
    }
    if (H5Tequal(type_id, H5T_NATIVE_INT) || H5Tequal(type_id, H5T_STD_I32LE))
    {
        if (H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids) == HDF5_error)
            *ierr = 5;
    }
    else if (H5Tequal(type_id, H5T_NATIVE_LLONG) || H5Tequal(type_id, H5T_STD_I64LE))
    {
        long long *t64 = (long long *)malloc((size_t)nlen * sizeof(long long));
        if (!t64)
            *ierr = 5;
        else
        {
            if (H5Dread(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, t64) == HDF5_error)
                *ierr = 5;
            else
                for (i = 0; i < nlen; i++)
                    ids[i] = (int)t64[i];
            free(t64);
        }
    }
    else
        *ierr = 5;
    H5Tclose(type_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Gclose(particles_id);
    H5Fclose(file_id);
}

/*
 * read required chemistry columns from sidecar into dat via Fortran callbacks
 */
void phantom_hdf5_extra_read(char *filename, int *ntotal, int *icomp_col_start, int *ncomp,
                             int *isrequired, int *ierr)
{
    hid_t file_id, particles_id;
    int nc = 0, np = 0;

    (void)ntotal;
    *ierr = 0;
    if (*ncomp <= 0)
        return;
    extra_icomp_col_start = *icomp_col_start;
    reset_extra_nonfinite_tracking(0);
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == HDF5_error)
    {
        *ierr = 1;
        return;
    }
    particles_id = open_group(file_id, PHANTOM_PARTICLES);
    if (particles_id == HDF5_error)
    {
        H5Fclose(file_id);
        *ierr = 2;
        return;
    }
    process_particles_group(particles_id, 1, 0, &nc, &np, isrequired);
    H5Gclose(particles_id);
    H5Fclose(file_id);
    print_extra_nonfinite_warning();
}
