#ifndef Model_OBJ_H_
#define Model_OBJ_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "Octree.h"
#include "BVH.h"
#include <map>
#include <cstdlib>
#include <igl/readOBJ.h>
using namespace std;
/*************************************************************************** 
  OBJ Loading 
 ***************************************************************************/
 
class Model_OBJ
{
  public: 
  	struct Edge_Info
  	{
  		int face_x, face_y;
  		map<int, int> loop_index;
  	};
  vector<set<int> > v_faces;

	Model_OBJ();			
  int Load(char *filename);	// Loads the model

 	void Calc_Bounding_Box();

  void Process_Manifold(int resolution);
  void Build_Tree(int resolution);
  void Build_BVH();
  void Construct_Manifold();
  void Project_Manifold();
  bool Project(glm::dvec3& o, glm::dvec3& d);
  void Save(const char* filename, bool color);
  void SaveOBJ(const char* filename);
  glm::dvec3 Closest_Point( const glm::dvec3 *triangle, const glm::dvec3 &sourcePosition );
  glm::dvec3 Find_Closest(int i);
  int is_manifold();
  bool Split_Grid(map<Grid_Index,int>& vcolor, vector<glm::dvec3>& nvertices, vector<glm::ivec4>& nface_indices, vector<set<int> >& v_faces, vector<glm::ivec3>& triangles);
  double clamp(double d1, double l, double r)
  {
    if (d1 < l)
      return l;
    if (d1 > r)
      return l;
    return d1;
  }
  vector<Grid_Index > v_info;
  
	vector<glm::dvec3> vertices, vertices_buf;
  vector<glm::dvec3> colors;
	vector<glm::ivec3> face_indices, face_indices_buf;
  vector<glm::dvec3> face_normals;
	
	glm::dvec3 min_corner, max_corner;
  Octree* tree;
  BVH* bvh;
  vector<BV*> bvs;
  char* fn;
  // vector field
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

};

#endif
