#define TINYOBJLOADER_IMPLEMENTATION
#include "Mesh.h"
#include <tiny_obj_loader.h>
#include "mctable.h"


int f_cnt = 1;

Mesh::Mesh()
	: m_uv(false), m_normal(false)
{

}

Mesh Mesh::LoadMesh(const std::string &filename, const std::string &newmesh)
{

	//ref: https://github.com/syoyo/tinyobjloader/tree/v0.9.x
	std::string inputfile = filename;
	std::string outputfile = newmesh;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	Mesh points;
	std::string err;
	bool ret = tinyobj::LoadObj(shapes, materials, err, inputfile.c_str());
	FILE *fp;
	fp = fopen(outputfile.c_str(), "w");
	fprintf(fp, "vt 1.000000 1.000000\n");
	fclose(fp);

	if (!err.empty()) { // `err` may contain warning message.
		std::cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	if (!fp) {
		std::cout << "Fail to open file: " << outputfile.c_str() << std::endl;
	}

	for (int i = 0;i < 20;++i) {
		for (int j = 0; j < 20; ++j) {
			for (int k = 0; k < 20; ++k) {
				points.pos[i][j][k] = 0;
			}
		}
	}

	std::cout << "# of shapes    : " << shapes.size() << std::endl;
	std::cout << "# of materials : " << materials.size() << std::endl;

	for (size_t i = 0; i < shapes.size(); i++) {
		printf("shape[%ld].name = %s\n", i, shapes[i].name.c_str());
		printf("Size of shape[%ld].indices: %ld\n", i, shapes[i].mesh.indices.size());
		printf("Size of shape[%ld].texcoords: %ld\n", i, shapes[i].mesh.texcoords.size());
		printf("Size of shape[%ld].material_ids: %ld\n", i, shapes[i].mesh.material_ids.size());

		printf("shape[%ld].vertices: %ld\n", i, shapes[i].mesh.positions.size());
		assert((shapes[i].mesh.indices.size() % 4) == 0);
		for (size_t f = 0; f < shapes[i].mesh.indices.size() / 4; f++) {
		/*
			printf("  idx[%ld] = %d, %d, %d, %d. mat_id = %d\n", f,
				shapes[i].mesh.indices[4 * f + 0],
				shapes[i].mesh.indices[4 * f + 1],
				shapes[i].mesh.indices[4 * f + 2],
				shapes[i].mesh.indices[4 * f + 3],
				shapes[i].mesh.material_ids[f]);
		*/
		}

		//v: position
		assert((shapes[i].mesh.positions.size() % 3) == 0);
		for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
		//for (size_t v = 0; v < 20; v++) {	//test
			/*
			printf("  v[%ld] = (%.2f, %.2f, %.2f)\n", v,
				shapes[i].mesh.positions[3 * v + 0],
				shapes[i].mesh.positions[3 * v + 1],
				shapes[i].mesh.positions[3 * v + 2]);
				
			fprintf(fp, "v %.1f %.1f %.1f\n",
				shapes[i].mesh.positions[3 * v + 0],
				shapes[i].mesh.positions[3 * v + 1],
				shapes[i].mesh.positions[3 * v + 2]);
			*/
			int x_index = (int)(shapes[i].mesh.positions[3 * v + 0] * 10 + 10);
			int y_index = (int)(shapes[i].mesh.positions[3 * v + 1] * 10 + 10);
			int z_index = (int)(shapes[i].mesh.positions[3 * v + 2] * 10 + 10);

			//printf("  v[%ld] = (%d, %d, %d)\n", v, x_index, y_index, z_index);
			points.pos[x_index][y_index][z_index] = 1;
			/*
			fprintf(fp, "v %.2f %.2f %.2f\n",
				shapes[i].mesh.positions[3 * v + 0],
				shapes[i].mesh.positions[3 * v + 1],
				shapes[i].mesh.positions[3 * v + 2]);
			*/
		}

		//vt: texcoords
		for (size_t v = 0; v < shapes[i].mesh.texcoords.size() / 3; v++) {
			/*
			printf("  vt[%ld] = (%f, %f)\n", v,
				shapes[i].mesh.texcoords[2 * v + 0],
				shapes[i].mesh.texcoords[2 * v + 1]);
			
			fprintf(fp, "vt %.2f %.2f \n",
				shapes[i].mesh.texcoords[2 * v + 0],
				shapes[i].mesh.texcoords[2 * v + 1]);
			*/
		}

		//vn: normals
		for (size_t v = 0; v < shapes[i].mesh.normals.size() / 3; v++) {
			/*
			printf("  vn[%ld] = (%f, %f, %f)\n", v,
				shapes[i].mesh.normals[3 * v + 0],
				shapes[i].mesh.normals[3 * v + 1],
				shapes[i].mesh.normals[3 * v + 2]);
			
			fprintf(fp, "vn %.2f %.2f %.2f\n",
				shapes[i].mesh.normals[3 * v + 0],
				shapes[i].mesh.normals[3 * v + 1],
				shapes[i].mesh.normals[3 * v + 2]);
			*/
		}

	}
	

	//material
	for (size_t i = 0; i < materials.size(); i++) {
		printf("material[%ld].name = %s\n", i, materials[i].name.c_str());
		printf("  material.Ka = (%f, %f ,%f)\n", materials[i].ambient[0], materials[i].ambient[1], materials[i].ambient[2]);
		printf("  material.Kd = (%f, %f ,%f)\n", materials[i].diffuse[0], materials[i].diffuse[1], materials[i].diffuse[2]);
		printf("  material.Ks = (%f, %f ,%f)\n", materials[i].specular[0], materials[i].specular[1], materials[i].specular[2]);
		printf("  material.Tr = (%f, %f ,%f)\n", materials[i].transmittance[0], materials[i].transmittance[1], materials[i].transmittance[2]);
		printf("  material.Ke = (%f, %f ,%f)\n", materials[i].emission[0], materials[i].emission[1], materials[i].emission[2]);
		printf("  material.Ns = %f\n", materials[i].shininess);
		printf("  material.Ni = %f\n", materials[i].ior);
		printf("  material.dissolve = %f\n", materials[i].dissolve);
		printf("  material.illum = %d\n", materials[i].illum);
		printf("  material.map_Ka = %s\n", materials[i].ambient_texname.c_str());
		printf("  material.map_Kd = %s\n", materials[i].diffuse_texname.c_str());
		printf("  material.map_Ks = %s\n", materials[i].specular_texname.c_str());
		printf("  material.map_Ns = %s\n", materials[i].specular_highlight_texname.c_str());
		std::map<std::string, std::string>::const_iterator it(materials[i].unknown_parameter.begin());
		std::map<std::string, std::string>::const_iterator itEnd(materials[i].unknown_parameter.end());
		for (; it != itEnd; it++) {
			printf("  material.%s = %s\n", it->first.c_str(), it->second.c_str());
		}
		printf("\n");

	}
	

	for (int i = 0;i < 19;++i) {
		for (int j = 0; j < 19; ++j) {
			for (int k = 0; k < 19; ++k) {

				GRIDCELL gridcell;
				TRIANGLE triangle[8];
				double isolevel = 0.5;

				for (int i = 0; i < 8; ++i) {
					gridcell.p[i].x = 0; gridcell.p[i].y = 0; gridcell.p[i].z = 0;
					gridcell.val[i]  = 0;
					triangle[i].p[0].x = 0; triangle[i].p[0].y = 0; triangle[i].p[0].z = 0;
					triangle[i].p[1].x = 0; triangle[i].p[1].y = 0; triangle[i].p[1].z = 0;
					triangle[i].p[2].x = 0; triangle[i].p[2].y = 0; triangle[i].p[2].z = 0;
				}
				
				gridcell.val[0] = points.pos[i][j][k];
				gridcell.p[0].x = (((double)(i)-10) / 10);
				gridcell.p[0].y = (((double)(j)-10) / 10);
				gridcell.p[0].z = (((double)(k)-10) / 10);
				gridcell.val[1] = points.pos[i+1][j][k];
				gridcell.p[1].x = (((double)(i+1)-10) / 10);
				gridcell.p[1].y = (((double)(j  )-10) / 10);
				gridcell.p[1].z = (((double)(k) - 10) / 10);
				gridcell.val[2] = points.pos[i + 1][j+1][k];
				gridcell.p[2].x = (((double)(i + 1) - 10) / 10);
				gridcell.p[2].y = (((double)(j+1)-10) / 10);
				gridcell.p[2].z = (((double)(k)-10) / 10);
				gridcell.val[3] = points.pos[i][j+1][k];
				gridcell.p[3].x = (((double)(i) - 10) / 10);
				gridcell.p[3].y = (((double)(j+1)-10) / 10);
				gridcell.p[3].z = (((double)(k) - 10) / 10);
				gridcell.val[4] = points.pos[i][j][k+1];
				gridcell.p[4].x = (((double)(i)-10) / 10);
				gridcell.p[4].y = (((double)(j) - 10) / 10);
				gridcell.p[4].z = (((double)(k+1)-10) / 10);
				gridcell.val[5] = points.pos[i+1][j][k + 1];
				gridcell.p[5].x = (((double)(i+1)-10) / 10);
				gridcell.p[5].y = (((double)(j) - 10) / 10);
				gridcell.p[5].z = (((double)(k + 1) - 10) / 10);
				gridcell.val[6] = points.pos[i + 1][j + 1][k+1];
				gridcell.p[6].x = (((double)(i + 1) - 10) / 10);
				gridcell.p[6].y = (((double)(j + 1) - 10) / 10);
				gridcell.p[6].z = (((double)(k+1)-10) / 10);
				gridcell.val[7] = points.pos[i][j + 1][k + 1];
				gridcell.p[7].x = (((double)(i) - 10) / 10);
				gridcell.p[7].y = (((double)(j + 1) - 10) / 10);
				gridcell.p[7].z = (((double)(k + 1) - 10) / 10);
				
				
				int ntriang = points.Polygonise(outputfile.c_str(), gridcell, isolevel, triangle);
				/*
				for(int i=0; i<ntriang; ++i){
					//printf("%lf %lf %lf\n", triangle[ntriang].p[0].x, triangle[ntriang].p[0].y, triangle[ntriang].p[0].z );
					//printf("%lf %lf %lf\n", triangle[ntriang].p[1].x, triangle[ntriang].p[1].y, triangle[ntriang].p[1].z );
					//printf("%lf %lf %lf\n", triangle[ntriang].p[2].x, triangle[ntriang].p[2].y, triangle[ntriang].p[2].z );
				}
				
				if(ntriang != 0)
					printf("ntriang:%d\n", ntriang);
				for (int i = 0; i < 8; ++i) 
					printf("gridcell[%d] val = %lf  \n", i, gridcell.val[i]);
				*/
				
			}
		}
	}
	
}

int Mesh::Polygonise(const std::string &filename, GRIDCELL grid, double isolevel, TRIANGLE  *triangle) {
	int ntriang = 0;
	int cubeindex = 0;
	XYZ vertlist[12];
	std::string outputfile = filename;
	FILE *fp;
	fp = fopen(outputfile.c_str(), "a");

	for (int i = 0;i < 12;++i) {
		vertlist[i].x = 0;
		vertlist[i].y = 0;
		vertlist[i].z = 0;
	}

	if (grid.val[0] > 0)	cubeindex |=   1;
	if (grid.val[1] > 0)	cubeindex |=   2;
	if (grid.val[2] > 0)	cubeindex |=   4;
	if (grid.val[3] > 0)	cubeindex |=   8;
	if (grid.val[4] > 0)	cubeindex |=  16;
	if (grid.val[5] > 0)	cubeindex |=  32;
	if (grid.val[6] > 0)	cubeindex |=  64;
	if (grid.val[7] > 0)	cubeindex |= 128;

	/* Cube is entirely in/out of the surface */
	int check_0 = 1;
	for (int i = 0;i < 6;++i) {
		if(mctable[cubeindex][i] != 0 )
			check_0 = 0;
	}
	if(check_0 == 1)
		return 0;
		
	/* Find the vertices where the surface intersects the cube */
	if (edgeTable[cubeindex] & 1)
		vertlist[0] =
		VertexInterp(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1]);
	if (edgeTable[cubeindex] & 2)
		vertlist[1] =
		VertexInterp(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
	if (edgeTable[cubeindex] & 4)
		vertlist[2] =
		VertexInterp(isolevel, grid.p[2], grid.p[3], grid.val[2], grid.val[3]);
	if (edgeTable[cubeindex] & 8)
		vertlist[3] =
		VertexInterp(isolevel, grid.p[3], grid.p[0], grid.val[3], grid.val[0]);
	if (edgeTable[cubeindex] & 16)
		vertlist[4] =
		VertexInterp(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
	if (edgeTable[cubeindex] & 32)
		vertlist[5] =
		VertexInterp(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
	if (edgeTable[cubeindex] & 64)
		vertlist[6] =
		VertexInterp(isolevel, grid.p[6], grid.p[7], grid.val[6], grid.val[7]);
	if (edgeTable[cubeindex] & 128)
		vertlist[7] =
		VertexInterp(isolevel, grid.p[7], grid.p[4], grid.val[7], grid.val[4]);
	if (edgeTable[cubeindex] & 256)
		vertlist[8] =
		VertexInterp(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
	if (edgeTable[cubeindex] & 512)
		vertlist[9] =
		VertexInterp(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
	if (edgeTable[cubeindex] & 1024)
		vertlist[10] =
		VertexInterp(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
	if (edgeTable[cubeindex] & 2048)
		vertlist[11] =
		VertexInterp(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7]);

	//create the triangle
	for (int i = 0;mctable[cubeindex][i] != -1;i += 3) {
		triangle[ntriang].p[0] = vertlist[mctable[cubeindex][i  ]];
		triangle[ntriang].p[1] = vertlist[mctable[cubeindex][i+1]];
		triangle[ntriang].p[2] = vertlist[mctable[cubeindex][i+2]];

		printf("%lf %lf %lf\n", triangle[ntriang].p[0].x, triangle[ntriang].p[0].y, triangle[ntriang].p[0].z );
		printf("%lf %lf %lf\n", triangle[ntriang].p[1].x, triangle[ntriang].p[1].y, triangle[ntriang].p[1].z );
		printf("%lf %lf %lf\n", triangle[ntriang].p[2].x, triangle[ntriang].p[2].y, triangle[ntriang].p[2].z );
		//v
		fprintf(fp, "v %lf %lf %lf\n", triangle[ntriang].p[0].x, triangle[ntriang].p[0].y, triangle[ntriang].p[0].z);
		fprintf(fp, "v %lf %lf %lf\n", triangle[ntriang].p[1].x, triangle[ntriang].p[1].y, triangle[ntriang].p[1].z);
		fprintf(fp, "v %lf %lf %lf\n", triangle[ntriang].p[2].x, triangle[ntriang].p[2].y, triangle[ntriang].p[2].z);
		
		//vn
		//resolution: 20
		//p[2]
		TRIANGLE vertex_normal = Vertex_normal(triangle[ntriang], ntriang);
		//printf("%lf %lf %lf\n", vertex_normal.p[0].x, vertex_normal.p[0].y, vertex_normal.p[0].z);
		//printf("%lf %lf %lf\n", vertex_normal.p[1].x, vertex_normal.p[1].y, vertex_normal.p[1].z);
		//printf("%lf %lf %lf\n", vertex_normal.p[2].x, vertex_normal.p[2].y, vertex_normal.p[2].z);
		fprintf(fp, "vn %lf %lf %lf\n", vertex_normal.p[0].x, vertex_normal.p[0].y, vertex_normal.p[0].z);
		fprintf(fp, "vn %lf %lf %lf\n", vertex_normal.p[1].x, vertex_normal.p[1].y, vertex_normal.p[1].z);
		fprintf(fp, "vn %lf %lf %lf\n", vertex_normal.p[2].x, vertex_normal.p[2].y, vertex_normal.p[2].z);


		//f
		fprintf(fp, "f %d/1/%d %d/1/%d %d/1/%d\n", f_cnt, f_cnt, f_cnt+1, f_cnt+1, f_cnt+2, f_cnt+2);
		f_cnt += 3;
		ntriang++;

	}

	fclose(fp);
	return ntriang;

}

TRIANGLE Mesh::Vertex_normal(TRIANGLE  triangle, int ntriang) {
	TRIANGLE vertex_normal;
	double vertex_normal_length = 0;
	//p[0]
	vertex_normal.p[0].x = 0.5 * (calcNormal(triangle.p[0].x + 1, triangle.p[0].y, triangle.p[0].z) \
								- calcNormal(triangle.p[0].x - 1, triangle.p[0].y, triangle.p[0].z)) / 20;
	vertex_normal.p[0].y = 0.5 * (calcNormal(triangle.p[0].x, triangle.p[0].y + 1, triangle.p[0].z) \
								- calcNormal(triangle.p[0].x, triangle.p[0].y - 1, triangle.p[0].z)) / 20;
	vertex_normal.p[0].z = 0.5 * (calcNormal(triangle.p[0].x, triangle.p[0].y, triangle.p[0].z + 1) \
								- calcNormal(triangle.p[0].x, triangle.p[0].y, triangle.p[0].z - 1)) / 20;
	vertex_normal_length = sqrt((vertex_normal.p[0].x * vertex_normal.p[0].x) + \
										(vertex_normal.p[0].y * vertex_normal.p[0].y) + \
										(vertex_normal.p[0].z * vertex_normal.p[0].z));
	if (vertex_normal_length != 0) {
		vertex_normal.p[0].x = vertex_normal.p[0].x / vertex_normal_length;
		vertex_normal.p[0].y = vertex_normal.p[0].y / vertex_normal_length;
		vertex_normal.p[0].z = vertex_normal.p[0].z / vertex_normal_length;
	}

	//p[1]
	vertex_normal.p[1].x = 0.5 * (calcNormal(triangle.p[1].x + 1, triangle.p[1].y, triangle.p[1].z) \
								- calcNormal(triangle.p[1].x - 1, triangle.p[1].y, triangle.p[1].z)) / 20;
	vertex_normal.p[1].y = 0.5 * (calcNormal(triangle.p[1].x, triangle.p[1].y + 1, triangle.p[1].z) \
								- calcNormal(triangle.p[1].x, triangle.p[1].y - 1, triangle.p[1].z)) / 20;
	vertex_normal.p[1].z = 0.5 * (calcNormal(triangle.p[1].x, triangle.p[1].y, triangle.p[1].z + 1) \
								- calcNormal(triangle.p[1].x, triangle.p[1].y, triangle.p[1].z - 1)) / 20;
	vertex_normal_length = sqrt((vertex_normal.p[1].x * vertex_normal.p[1].x) + \
										(vertex_normal.p[1].y * vertex_normal.p[1].y) + \
										(vertex_normal.p[1].z * vertex_normal.p[1].z));
	if (vertex_normal_length != 0) {
		vertex_normal.p[1].x = vertex_normal.p[1].x / vertex_normal_length;
		vertex_normal.p[1].y = vertex_normal.p[1].y / vertex_normal_length;
		vertex_normal.p[1].z = vertex_normal.p[1].z / vertex_normal_length;
	}

	//p[2]
	vertex_normal.p[2].x = 0.5 * (calcNormal(triangle.p[2].x + 1, triangle.p[2].y, triangle.p[2].z) \
								- calcNormal(triangle.p[2].x - 1, triangle.p[2].y, triangle.p[2].z)) / 20;
	vertex_normal.p[2].y = 0.5 * (calcNormal(triangle.p[2].x, triangle.p[2].y + 1, triangle.p[2].z) \
								- calcNormal(triangle.p[2].x, triangle.p[2].y - 1, triangle.p[2].z)) / 20;
	vertex_normal.p[2].z = 0.5 * (calcNormal(triangle.p[2].x, triangle.p[2].y, triangle.p[2].z + 1) \
								- calcNormal(triangle.p[2].x, triangle.p[2].y, triangle.p[2].z - 1)) / 20;
	vertex_normal_length = sqrt((vertex_normal.p[2].x * vertex_normal.p[2].x) + \
										(vertex_normal.p[2].y * vertex_normal.p[2].y) + \
										(vertex_normal.p[2].z * vertex_normal.p[2].z));
	if (vertex_normal_length != 0) {
		vertex_normal.p[2].x = vertex_normal.p[2].x / vertex_normal_length;
		vertex_normal.p[2].y = vertex_normal.p[2].y / vertex_normal_length;
		vertex_normal.p[2].z = vertex_normal.p[2].z / vertex_normal_length;
	}
	return vertex_normal;
}

XYZ Mesh::VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2) {
	double mu;
	XYZ p;

	if (fabs(isolevel - valp1) < 0.00001)
		return(p1);
	if (fabs(isolevel - valp2) < 0.00001)
		return(p2);
	if (fabs(valp1 - valp2) < 0.00001)
		return(p1);
	mu = (isolevel - valp1) / (valp2 - valp1);
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);

	return(p);

}


//ref: https://github.com/lettier/isosurface/blob/master/source/assets/scripts/index.js
double Mesh::calcNormal(double x, double y, double z) {
	double x2 = x*x;
	double y2 = y*y;
	double z2 = z*z;

	double x4 = x2*x2;
	double y4 = y2*y2;
	double z4 = z2*z2;

	double a = -1.0;
	double b =  0.0;
	double c =  0.5;

	double d = x2 + y2 + z2;
	double d2 = d*d;

	double value = x4 + y4 + z4 + a * d2 + b * d + c;

	return value;
}

void Mesh::release()
{
	glDeleteVertexArrays(1, &vao);
	glDeleteBuffers(3, vbo);
	glDeleteBuffers(1, &ibo);

}

void Mesh::draw()
{
	glBindVertexArray(vao);
	glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, nullptr);
}

bool Mesh::hasNormal() const
{
	return m_normal;
}

bool Mesh::hasUV() const
{
	return m_uv;
}
