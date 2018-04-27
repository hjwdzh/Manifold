#include "igl/circulation.h"
#include "igl/collapse_edge.h"
#include "igl/edge_flaps.h"
#include "igl/shortest_edge_and_midpoint.h"
#include "igl/read_triangle_mesh.h"
#include "igl/write_triangle_mesh.h"
#include "igl/decimate.h"
#include "igl/qslim.h"
#include <Eigen/Core>
#include <iostream>

int main(int argc, char * argv[])
{
  if (argc == 1) {
    printf("./simplify -i input.obj -o output.obj -m -c 1e-2 -f 10000 -r 0.2\n");
    return 0;
  }
  using namespace std;
  using namespace Eigen;
  using namespace igl;

  string input_name, output_name;
  double max_cost = 1e30;
  double max_ratio = 2.;
  int max_faces = 0x7fffffff;
  bool check_manifold = false;

  for (int i = 0; i < argc; ++i) {
    if (i + 1 < argc) {
      if (strcmp(argv[i], "-i") == 0) {
        input_name = argv[i + 1];
      }
      if (strcmp(argv[i], "-o") == 0) {
        output_name = argv[i + 1];
      }
      if (strcmp(argv[i], "-f") == 0) {
        sscanf(argv[i + 1], "%d", &max_faces);
      }
      if (strcmp(argv[i], "-r") == 0) {
        sscanf(argv[i + 1], "%lf", &max_ratio);
      }
      if (strcmp(argv[i], "-c") == 0) {
        sscanf(argv[i + 1], "%lf", &max_cost);
      }
    }
    if (strcmp(argv[i], "-m") == 0) {
      check_manifold = true;
    }
  }

  MatrixXd V,OV;
  MatrixXi F,OF;
  VectorXi J, I;
  read_triangle_mesh(input_name,OV,OF);

  if (max_faces == 0x7fffffff) {
    max_faces = max_ratio * OF.rows() + 0.5;
  } else {
    if (max_ratio < 1)
      max_faces = std::max(max_faces, (int)(max_ratio * OF.rows() + 0.5));
  }

  auto MyDecimate = [&](
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const size_t max_m,
    const double max_cost,
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & G,
    Eigen::VectorXi & J,
    Eigen::VectorXi & I) {
    // Original number of faces
    const int orig_m = F.rows();
    // Tracking number of faces
    int m = F.rows();
    typedef Eigen::MatrixXd DerivedV;
    typedef Eigen::MatrixXi DerivedF;
    DerivedV VO;
    DerivedF FO;
    igl::connect_boundary_to_infinity(V,F,VO,FO);
    // decimate will not work correctly on non-edge-manifold meshes. By extension
    // this includes meshes with non-manifold vertices on the boundary since these
    // will create a non-manifold edge when connected to infinity.
    if(!is_edge_manifold(FO))
    {
      return false;
    }

    auto stopping_condition = 
    [&](
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &E,
    const Eigen::VectorXi &EMAP,
    const Eigen::MatrixXi &EF,
    const Eigen::MatrixXi &EI,
    const std::set<std::pair<double,int> > &Q,
    const std::vector<std::set<std::pair<double,int> >::iterator > &QIT,
    const Eigen::MatrixXd &C,
    const int e,
    const int e1,
    const int e2,
    const int f1,
    const int f2)->bool
    {
      // Only subtract if we're collapsing a real face
      if(f1 < orig_m) m-=1;
      if(f2 < orig_m) m-=1;
      return (m<=(int)max_m) || (Q.begin()->first >= max_cost);
    };

    using namespace igl;

    Eigen::VectorXi EMAP;
    Eigen::MatrixXi E,EF,EI;
    edge_flaps(FO,E,EMAP,EF,EI);
    // Quadrics per vertex
    typedef std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> Quadric;
    std::vector<Quadric> quadrics;
    per_vertex_point_to_plane_quadrics(VO,FO,EMAP,EF,EI,quadrics);
    // State variables keeping track of edge we just collapsed
    int v1 = -1;
    int v2 = -1;
    // Callbacks for computing and updating metric
    std::function<void(
      const int e,
      const Eigen::MatrixXd &,
      const Eigen::MatrixXi &,
      const Eigen::MatrixXi &,
      const Eigen::VectorXi &,
      const Eigen::MatrixXi &,
      const Eigen::MatrixXi &,
      double &,
      Eigen::RowVectorXd &)> cost_and_placement;
    std::function<bool(
      const Eigen::MatrixXd &                                         ,/*V*/
      const Eigen::MatrixXi &                                         ,/*F*/
      const Eigen::MatrixXi &                                         ,/*E*/
      const Eigen::VectorXi &                                         ,/*EMAP*/
      const Eigen::MatrixXi &                                         ,/*EF*/
      const Eigen::MatrixXi &                                         ,/*EI*/
      const std::set<std::pair<double,int> > &                        ,/*Q*/
      const std::vector<std::set<std::pair<double,int> >::iterator > &,/*Qit*/
      const Eigen::MatrixXd &                                         ,/*C*/
      const int                                                        /*e*/
      )> pre_collapse;
    std::function<void(
      const Eigen::MatrixXd &                                         ,   /*V*/
      const Eigen::MatrixXi &                                         ,   /*F*/
      const Eigen::MatrixXi &                                         ,   /*E*/
      const Eigen::VectorXi &                                         ,/*EMAP*/
      const Eigen::MatrixXi &                                         ,  /*EF*/
      const Eigen::MatrixXi &                                         ,  /*EI*/
      const std::set<std::pair<double,int> > &                        ,   /*Q*/
      const std::vector<std::set<std::pair<double,int> >::iterator > &, /*Qit*/
      const Eigen::MatrixXd &                                         ,   /*C*/
      const int                                                       ,   /*e*/
      const int                                                       ,  /*e1*/
      const int                                                       ,  /*e2*/
      const int                                                       ,  /*f1*/
      const int                                                       ,  /*f2*/
      const bool                                                  /*collapsed*/
      )> post_collapse;
    qslim_optimal_collapse_edge_callbacks(
      E,quadrics,v1,v2, cost_and_placement, pre_collapse,post_collapse);
    // Call to greedy decimator
    bool ret = decimate(
      VO, FO,
      cost_and_placement,
      stopping_condition,
      pre_collapse,
      post_collapse,
      E, EMAP, EF, EI,
      U, G, J, I);
    // Remove phony boundary faces and clean up
    const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
    igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
    igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
    Eigen::VectorXi _1,I2;
    igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
    igl::slice(Eigen::VectorXi(I),I2,1,I);

    return ret;
  };

  MyDecimate(OV, OF, max_faces, max_cost, V, F, J, I);

  if (check_manifold) {
    std::set<std::pair<int, int> > dedges;
    for (int i = 0; i < F.rows(); ++i) {
      for (int j = 0; j < 3; ++j) {
        int v1 = F(i, j);
        int v2 = F(i, (j + 1) % 3);
        auto key = std::make_pair(v1, v2);
        if (dedges.count(key)) {
          printf("Non manifold!\n");
          return 0;
        }
        dedges.insert(key);
      }
    }
    printf("Output is a manifold!\n");
  }

  write_triangle_mesh(output_name,V,F);

  return 0;
}
