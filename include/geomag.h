
#ifndef GEOMAG_H
#define GEOMAG_H

#ifdef __cplusplus
extern "C" {
#endif

/* Max number of models in a file */
#define MAXMOD 30

typedef enum {
	kUnitsKilometers=1,
	kUnitsMeters    =2,
	kUnitsFeet      =3
} Units;

typedef enum {
	kCoordSysGeodetic  =1, /* WGS-84 */
	kCoordSysGeocentric=2
} CoordinateSystem;

typedef struct {
	double d;    /* Declination of the field from the
	                geographic north (deg). */
	double i;    /* Inclination of the field (deg). */
	double h;
	double f;    /* Total field intensity. */
	double x;
	double y;
	double z;
	double ddot; /* Annual rate of change of declination. (arc-min/yr) */
	double fdot; /* Annual rate of change of . */
	double hdot; /* Annual rate of change of . */
	double idot; /* Annual rate of change of inclination. (arc-min/yr). */
	double xdot; /* Annual rate of change of . */
	double ydot; /* Annual rate of change of . */
	double zdot; /* Annual rate of change of . */
} BField;

#define MAXDEG 13
#define MAXCOEFF (MAXDEG*(MAXDEG+2)+1) /* index starts with 1!, (from old Fortran?) */

typedef struct {
	char name[MAXMOD][9];
	double epoch[MAXMOD];  /* epoch of model. */
	double yrmin[MAXMOD];  /* Min year of model. */
	double yrmax[MAXMOD];  /* Max year of model. */
	double minyr;          /* Min year of all models. */
	double maxyr;          /* Max year of all models. */
	double altmin[MAXMOD]; /* Minimum height of each model. */
	double altmax[MAXMOD]; /* Maximum height of each model. */
	int nmodel;            /* Number of models in file */
	int max1[MAXMOD];      /* Main field coefficient. */
	int max2[MAXMOD];      /* Secular variation coefficient. */
	int max3[MAXMOD];      /* Acceleration coefficient. */
	long irec_pos[MAXMOD]; /* Record counter for header */
	double gh1[MAXMOD][MAXCOEFF]; /* Schmidt quasi-normal internal spherical harmonic coeff. */
	double gh2[MAXMOD][MAXCOEFF]; /* Schmidt quasi-normal internal spherical harmonic coeff. */
} BFieldModel;

extern int get_field_components(BField *const bfield,
                                BFieldModel const*const model,
                                double alt,
                                const Units altUnits,
                                const CoordinateSystem coordSys,
                                const double latitude,
                                const double longitude,
                                const double sdate);

extern double julday(const int month, const int day, const int year);

#ifndef TARGET_EMBEDDED
extern int read_model(BFieldModel *const model, const char mdfile[]);
#endif

extern BFieldModel mag_model; /* GeoMag model for compute declination */
extern BField      mag;       /* geomagnetic field values */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
