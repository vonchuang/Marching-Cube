#ifndef _STATIC_MESH_H_
#define _STATIC_MESH_H_

#include <GL/glew.h>
#include <string>
#include <iostream>
#include <ostream>
class StaticMesh {
public:
    void release();

    static StaticMesh LoadMesh(const std::string &filename);
    void draw();

	bool hasNormal() const;
	bool hasUV() const;
private:
    StaticMesh();
    GLuint vao;
    GLuint vbo[3];
    GLuint ibo;
    GLuint numIndices;

	bool m_uv, m_normal;
};

#endif