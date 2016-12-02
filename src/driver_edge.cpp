#include <apfBox.h>
#include <apf.h>
#include <PCU.h>
#include <iostream>

int main()
{
  MPI_Init(0, 0);
  PCU_Comm_Init();
  auto mesh = apf::makeMdsBox( 2, 2, 2, 1, 1, 1, true);

  apf::MeshEntity* edge;
  auto edges = mesh->begin(1);
    while ((edge = mesh->iterate(edges))) {
            apf::Downward vertices;
            int numVertices;
            numVertices = mesh->getDownward(edge, 0, vertices);
                  for (int i = 0; i < numVertices; ++i) {
                       apf::Vector3 point;
                       mesh->getPoint(vertices[i], 0, point);
                       std::cout << point << std::endl;
                  }
    }

  apf::writeVtkFiles( "Meshed_Box", mesh);
  PCU_Comm_Free();
  MPI_Finalize();

}
