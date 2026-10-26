/*
 * C implementation of the falcON HDF5 snapshot reader.
 *
 * Replaces the original C++ version (read_data_falcON_hdf5_utils.cc)
 * so that the reader can be linked into the main splash binary using
 * only the HDF5 C library (no -lhdf5_cpp or C++ compiler required).
 *
 * The extern "C" interface is identical to the original, so the
 * Fortran wrapper (read_data_falcON_hdf5.f90) needs no changes.
 *
 * Copyright (C) 2015  Walter Dehnen (original C++ version)
 * Copyright (C) 2026  Daniel Price  (C rewrite)
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <hdf5.h>

/*-----------------------------------------------------------------------
 *  Forward declarations of public functions (avoid implicit-declaration
 *  errors when one public function calls another before its definition)
 *-----------------------------------------------------------------------*/
void close_falcON_file(void);

/*-----------------------------------------------------------------------
 *  Constants
 *-----------------------------------------------------------------------*/
#define MAX_NUM_TYPES   6
#define MAX_FIELDS    128
#define MAX_NAME       64
#define MAX_PARTS       3   /* sink, gas, std */
#define PREFERRED_N    16

/*-----------------------------------------------------------------------
 *  Fortran callbacks (implemented in read_data_falcON_hdf5.f90)
 *-----------------------------------------------------------------------*/
void read_falcON_data_into_splash(const int *icol, const int *ndat,
                                  const double *data, const int *type);
void set_splash_block_label(const int *icol, const char *name);
void set_splash_particle_label(const int *ityp, const char *name);

/*-----------------------------------------------------------------------
 *  Per-particle-type bookkeeping
 *-----------------------------------------------------------------------*/
typedef struct {
    hid_t   group_id;       /* HDF5 group handle (>= 0 when open) */
    char    name[MAX_NAME]; /* "sink", "gas", or "std"             */
    int     number;         /* particle count                      */
} ParticleSet;

/*-----------------------------------------------------------------------
 *  Per-field bookkeeping (field name -> dimensionality)
 *-----------------------------------------------------------------------*/
typedef struct {
    char    name[MAX_NAME];
    hsize_t dims;           /* 1 = scalar, 3 = vector, 6 = tensor */
} FieldEntry;

/*-----------------------------------------------------------------------
 *  Module-level state
 *-----------------------------------------------------------------------*/
static int         Debug      = 0;
static hid_t       FileId     = -1;       /* currently open HDF5 file   */
static int         IndexSnap  = 0;        /* next snapshot index        */
static int         NumSnap    = 0;        /* snapshots in file          */
static int         SplashCol  = 0;        /* current splash column      */
static int         NumCol     = 0;        /* total splash columns       */
static ParticleSet Parts[MAX_PARTS];
static int         NumParts   = 0;
static FieldEntry  Fields[MAX_FIELDS];
static int         NumFields  = 0;
static int         falcON_read_error = 0;

/*-----------------------------------------------------------------------
 *  Helpers
 *-----------------------------------------------------------------------*/

/* look up a field by name; return index or -1 */
static int find_field(const char *name)
{
    int i;
    for (i = 0; i < NumFields; i++) {
        if (strcmp(Fields[i].name, name) == 0) return i;
    }
    return -1;
}

/* add a new field; return index */
static int add_field(const char *name, hsize_t dims)
{
    if (NumFields >= MAX_FIELDS) return -1;
    strncpy(Fields[NumFields].name, name, MAX_NAME - 1);
    Fields[NumFields].name[MAX_NAME - 1] = '\0';
    Fields[NumFields].dims = dims;
    return NumFields++;
}

/* close the currently open snapshot (groups + state) */
static void close_snapshot(void)
{
    int i;
    for (i = 0; i < NumParts; i++) {
        if (Parts[i].group_id >= 0) {
            H5Gclose(Parts[i].group_id);
            Parts[i].group_id = -1;
        }
    }
    NumParts  = 0;
    NumFields = 0;
    SplashCol = 0;
    NumCol    = 0;
}

/* field selection: should this field be read for this particle type? */
static int select_field(const char *field, const char *ptype)
{
    if (strcmp(field, "spin") == 0 || strcmp(field, "eabs") == 0 ||
        strcmp(field, "maxA") == 0)
        return (strcmp(ptype, "sink") == 0);

    if (strcmp(field, "snum") == 0 || strcmp(field, "uin")  == 0 ||
        strcmp(field, "entr") == 0 || strcmp(field, "dlKt") == 0 ||
        strcmp(field, "dlKe") == 0 || strcmp(field, "srho") == 0 ||
        strcmp(field, "alfa") == 0 || strcmp(field, "divv") == 0 ||
        strcmp(field, "dlht") == 0 || strcmp(field, "vsig") == 0 ||
        strcmp(field, "fact") == 0 || strcmp(field, "csnd") == 0 ||
        strcmp(field, "pres") == 0 || strcmp(field, "vort") == 0 ||
        strcmp(field, "dtdv") == 0 || strcmp(field, "qmin") == 0 ||
        strcmp(field, "delE") == 0 || strcmp(field, "coll") == 0)
        return (strcmp(ptype, "gas") == 0);

    if (strcmp(field, "krnH") == 0 || strcmp(field, "maxR") == 0)
        return (strcmp(ptype, "sink") == 0 || strcmp(ptype, "gas") == 0);

    return 1; /* everything else: read for all types */
}

/* finish one splash column: set its label and advance counter */
static void close_column(const char *name)
{
    if (SplashCol >= NumCol) {
        if (Debug) fprintf(stderr, "falcON: column overflow (%d >= %d)\n",
                           SplashCol, NumCol);
        return;
    }
    set_splash_block_label(&SplashCol, name);
    SplashCol++;
}

/*-----------------------------------------------------------------------
 *  Read one component of a field into the current splash column
 *-----------------------------------------------------------------------*/
static void read_column(const char *field, hsize_t comp, hsize_t dims)
{
    int     type, read_any = 0;
    char    colname[MAX_NAME];
    double *buf = NULL;

    for (type = 0; type < NumParts; type++) {
        hid_t   dset, space, memspace;
        hsize_t count[2], start[2];
        int     ndat, rank;

        if (!select_field(field, Parts[type].name)) continue;

        /* try to open the dataset; skip silently if absent */
        dset = H5Dopen2(Parts[type].group_id, field, H5P_DEFAULT);
        if (dset < 0) continue;

        space = H5Dget_space(dset);
        rank  = H5Sget_simple_extent_ndims(space);
        H5Sget_simple_extent_dims(space, count, NULL);
        ndat = (int)Parts[type].number;

        {
            double *newbuf = (double *)realloc(buf, (size_t)ndat * sizeof(double));
            if (!newbuf && ndat > 0) {
                free(buf);
                falcON_read_error = 1;
                H5Sclose(space);
                H5Dclose(dset);
                return;
            }
            buf = newbuf;
        }

        if (dims == 1) {
            /* scalar field: read the whole thing */
            H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, buf);
            /* halve kernel smoothing length for gas */
            if (strcmp(field, "krnH") == 0 && strcmp(Parts[type].name, "gas") == 0) {
                int k;
                for (k = 0; k < ndat; k++) buf[k] *= 0.5;
            }
        } else {
            /* vector/tensor component: hyperslab select */
            start[0] = 0;     start[1] = comp;
            count[0] = (hsize_t)ndat; count[1] = 1;
            H5Sselect_hyperslab(space, H5S_SELECT_SET, start, NULL,
                                count, NULL);
            memspace = H5Screate_simple(1, count, NULL);
            H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space,
                    H5P_DEFAULT, buf);
            H5Sclose(memspace);
        }

        if (Debug > 2)
            fprintf(stderr, "falcON: reading %d '%s' for '%s' (type=%d) "
                    "into col=%d\n", ndat, field, Parts[type].name,
                    type, SplashCol);

        read_falcON_data_into_splash(&SplashCol, &ndat, buf, &type);
        read_any = 1;

        /* krnH produces per-type column names */
        if (strcmp(field, "krnH") == 0) {
            if (strcmp(Parts[type].name, "gas") == 0) {
                close_column("h");
                read_any = 0;
            } else if (strcmp(Parts[type].name, "sink") == 0) {
                close_column("snkR");
                read_any = 0;
            }
        }

        H5Sclose(space);
        H5Dclose(dset);
    }

    free(buf);

    /* build column name for non-special fields */
    if (read_any) {
        if (strcmp(field, "pos") == 0) {
            colname[0] = "xyz"[comp];
            colname[1] = '\0';
        } else if (dims == 3) {
            colname[0] = field[0];
            colname[1] = "xyz"[comp];
            colname[2] = '\0';
        } else if (dims == 6) {
            const char *suffixes[] = {"xx","xy","xz","yy","yz","zz"};
            colname[0] = (char)toupper((unsigned char)field[0]);
            strcpy(colname + 1, suffixes[comp]);
        } else {
            strncpy(colname, field, MAX_NAME - 1);
            colname[MAX_NAME - 1] = '\0';
        }
        close_column(colname);
    }
}

/* read all components of a field */
static void read_field_entry(const FieldEntry *f)
{
    hsize_t comp;
    for (comp = 0; comp < f->dims; comp++)
        read_column(f->name, comp, f->dims);
}

/* try to read a field by name */
static void read_field_by_name(const char *name)
{
    int idx = find_field(name);
    if (idx < 0) {
        if (Debug)
            fprintf(stderr, "WARNING: falcON field '%s' not present\n", name);
        return;
    }
    read_field_entry(&Fields[idx]);
}

/*=======================================================================
 *  Public interface (called from Fortran via iso_c_binding)
 *=======================================================================*/

void set_falcON_debugging_level(const int *d)
{
    Debug = *d;
}

void open_falcON_file(const char *filename, int *ierr)
{
    hid_t attr;
    unsigned int ns;

    *ierr = 1;
    IndexSnap = 0;
    close_snapshot();

    if (Debug < 2) H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    if (filename == NULL || filename[0] == '\0') {
        if (Debug) fprintf(stderr, "open_falcON_file(): empty filename\n");
        close_falcON_file();
        return;
    }

    /* open HDF5 file */
    FileId = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (FileId < 0) {
        if (Debug)
            fprintf(stderr, "open_falcON_file(): cannot open '%s'\n", filename);
        return;
    }

    /* verify this is a falcON snapshot by checking for the attribute */
    attr = H5Aopen(FileId, "falcON", H5P_DEFAULT);
    if (attr < 0) {
        if (Debug)
            fprintf(stderr, "open_falcON_file(): '%s' is not a falcON file\n",
                    filename);
        H5Fclose(FileId);
        FileId = -1;
        return;
    }
    H5Aclose(attr);

    /* read number of snapshots */
    attr = H5Aopen(FileId, "num_snapshots", H5P_DEFAULT);
    if (attr < 0) {
        if (Debug)
            fprintf(stderr, "open_falcON_file(): missing num_snapshots\n");
        H5Fclose(FileId);
        FileId = -1;
        return;
    }
    if (H5Aread(attr, H5T_NATIVE_UINT32, &ns) < 0) {
        H5Aclose(attr);
        H5Fclose(FileId);
        FileId = -1;
        return;
    }
    H5Aclose(attr);
    NumSnap = (int)ns;
    *ierr = 0;
}

int falcON_file_is_open(void)
{
    return (FileId >= 0) ? 1 : 0;
}

/*
 * return 1 if filename is a falcon HDF5 snapshot (root attribute "falcON")
 */
int falcon_is_falcon_file(const char *filename)
{
    hid_t file_id, attr;
    int is_falcon = 0;

    if (filename == NULL || filename[0] == '\0') return 0;

    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) return 0;

    attr = H5Aopen(file_id, "falcON", H5P_DEFAULT);
    if (attr >= 0) {
        is_falcon = 1;
        H5Aclose(attr);
    }
    H5Fclose(file_id);
    return is_falcon;
}

void close_falcON_file(void)
{
    close_snapshot();
    if (FileId >= 0) {
        H5Fclose(FileId);
        FileId = -1;
    }
    IndexSnap = 0;
    NumSnap   = 0;
}

int have_falcON_snapshot(int *ierr)
{
    *ierr = (FileId < 0) ? 1 : 0;
    return IndexSnap < NumSnap;
}

int num_falcON_snapshots(int *ierr)
{
    *ierr = (FileId < 0) ? 1 : 0;
    return NumSnap;
}

/*-----------------------------------------------------------------------
 *  Callback used by H5Literate to count datasets in a particle group
 *-----------------------------------------------------------------------*/
struct iterate_data {
    hid_t  group_id;
    int   *ncol;
};

static herr_t count_fields_cb(hid_t group_id, const char *name,
                               const H5L_info_t *info, void *opdata)
{
    struct iterate_data *d = (struct iterate_data *)opdata;
    hid_t  dset, space;
    int    rank, idx;
    hsize_t dims_arr[2];
    hsize_t field_dims;

    (void)info; /* unused */

    dset = H5Dopen2(group_id, name, H5P_DEFAULT);
    if (dset < 0) return 0;

    space = H5Dget_space(dset);
    rank  = H5Sget_simple_extent_ndims(space);
    H5Sget_simple_extent_dims(space, dims_arr, NULL);
    field_dims = (rank == 1) ? 1 : dims_arr[1];
    H5Sclose(space);
    H5Dclose(dset);

    idx = find_field(name);
    if (idx >= 0) {
        /* krnH generates an extra column per particle type */
        if (strcmp(name, "krnH") == 0) (*(d->ncol))++;
    } else {
        add_field(name, field_dims);
        (*(d->ncol)) += (int)field_dims;
    }

    if (Debug > 2)
        fprintf(stderr, "falcON: counting field '%s' dims=%d ncol->%d\n",
                name, (int)field_dims, *(d->ncol));

    return 0;
}

void open_falcON_snapshot(int *ntype, int npart[MAX_NUM_TYPES],
                          int *ncol, int *dimX, int *dimV,
                          double *timeval, double hper[3], int *ierr)
{
    char snapname[32];
    hid_t snap_id, attr, part_grp;
    unsigned int number;
    int i;
    const char *types[MAX_PARTS] = {"sink", "gas", "std"};
    struct iterate_data idata;

    *ierr  = 1;
    *ntype = 0;
    *ncol  = 0;
    *dimX  = 3;
    *dimV  = 3;
    close_snapshot();

    if (FileId < 0 || IndexSnap >= NumSnap) return;

    sprintf(snapname, "snapshot%d", IndexSnap++);

    snap_id = H5Gopen2(FileId, snapname, H5P_DEFAULT);
    if (snap_id < 0) return;

    /* read time attribute */
    attr = H5Aopen(snap_id, "time", H5P_DEFAULT);
    if (attr < 0) {
        H5Gclose(snap_id);
        return;
    }
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, timeval) < 0) {
        H5Aclose(attr);
        H5Gclose(snap_id);
        return;
    }
    H5Aclose(attr);

    /* read hper attribute (periodic half-periods) */
    attr = H5Aopen(snap_id, "hper", H5P_DEFAULT);
    if (attr < 0) {
        H5Gclose(snap_id);
        return;
    }
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, hper) < 0) {
        H5Aclose(attr);
        H5Gclose(snap_id);
        return;
    }
    H5Aclose(attr);

    /* iterate over particle types */
    for (i = 0; i < MAX_PARTS; i++) {
        char nname[16];
        sprintf(nname, "N%s", types[i]);

        number = 0;
        attr = H5Aopen(snap_id, nname, H5P_DEFAULT);
        if (attr >= 0) {
            if (H5Aread(attr, H5T_NATIVE_UINT32, &number) < 0) number = 0;
            H5Aclose(attr);
        }
        if (number == 0) continue;

        part_grp = H5Gopen2(snap_id, types[i], H5P_DEFAULT);
        if (part_grp < 0) continue;

        Parts[NumParts].group_id = part_grp;
        strncpy(Parts[NumParts].name, types[i], MAX_NAME - 1);
        Parts[NumParts].name[MAX_NAME - 1] = '\0';
        Parts[NumParts].number = (int)number;
        npart[*ntype] = (int)number;
        (*ntype)++;
        NumParts++;

        /* count fields and columns */
        idata.group_id = part_grp;
        idata.ncol = ncol;
        H5Literate(part_grp, H5_INDEX_NAME, H5_ITER_NATIVE,
                   NULL, count_fields_cb, &idata);
    }

    H5Gclose(snap_id);
    NumCol = *ncol;
    *ierr  = 0;
}

void read_falcON_snapshot(int *ierr)
{
    /* preferred field order (same as original) */
    static const char *preferred[] = {
        "pos", "key", "vel", "acc", "mass", "pot", "pex", "rung",
        "krnH", "srho", "uin", "entr", "divv", "dlKt", "dlKe", "alfa"
    };
    int  i, j, npreferred;
    int *done; /* flags for fields already read */
    int  type;

    falcON_read_error = 0;
    *ierr = 0;
    npreferred = (int)(sizeof(preferred) / sizeof(preferred[0]));

    done = (int *)calloc((size_t)NumFields, sizeof(int));
    if (!done) {
        *ierr = 1;
        return;
    }

    /* read preferred fields first */
    for (i = 0; i < npreferred; i++) {
        j = find_field(preferred[i]);
        if (j >= 0 && !done[j]) {
            read_field_entry(&Fields[j]);
            done[j] = 1;
        }
    }

    /* read remaining fields */
    for (j = 0; j < NumFields; j++) {
        if (!done[j]) read_field_entry(&Fields[j]);
    }

    free(done);

    if (falcON_read_error) *ierr = 1;

    /* set particle type labels */
    for (type = 0; type < NumParts; type++)
        set_splash_particle_label(&type, Parts[type].name);
}
