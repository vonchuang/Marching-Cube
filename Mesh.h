#ifndef _MESH_H_
#define _MESH_H_

#include <GL/glew.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "struct.h"
class Mesh {
public:
	void release();
	int pos[20][20][20];
	static Mesh LoadMesh(const std::string &filename, const std::string &newmesh);
	void draw();
	double calcNormal(double x, double y, double z);
	int Polygonise(const std::string &filename, GRIDCELL grid,double isolevel, TRIANGLE  *triangle);
	XYZ VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2);
	TRIANGLE Vertex_normal(TRIANGLE  triangle, int ntriang);
	bool hasNormal() const;
	bool hasUV() const;
private:
	Mesh();
	GLuint vao;
	GLuint vbo[3];
	GLuint ibo;
	GLuint numIndices;

	
	bool m_uv, m_normal;
};

#endif