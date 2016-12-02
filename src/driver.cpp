#include <apfBox.h>
#include <apf.h>
#include <PCU.h>
#include <iostream>

int main()
{
  MPI_Init(0, 0);
  PCU_Comm_Init();
  auto mesh = apf::makeMdsBox( 4, 4, 4, 4, 4, 4, true);
  apf::MeshEntity* vertex; 
  auto vertices = mesh->begin(0);
  apf::Vector3 point;
  while ((vertex = mesh->iterate(vertices)))
  {
    mesh->getPoint(vertex, 0, point);
    std::cout << point << std::endl;
  }
  apf::writeVtkFiles( "Meshed_Box", mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}
