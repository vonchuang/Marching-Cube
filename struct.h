#pragma once
typedef struct _point3f {
	double x;
	double y;
	double z;
} XYZ;

typedef struct {
	XYZ p[3];
} TRIANGLE;

typedef struct {
	XYZ p[8];
	double val[8];
} GRIDCELL;



