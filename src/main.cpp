#include <stdlib.h>
#include <stdio.h>

#include "Model_OBJ.h"

int main(int argc, char** argv)
{
  Model_OBJ obj;
  int resolution = 20000;
  if (argc < 3)
  {
    cout << "./manifold input.obj output.obj [resolution=20000]\n";
    return 0;
  }
  obj.Load(argv[1]);

  if (argc > 3)
  {
    sscanf(argv[3], "%d", &resolution);
  }
  printf("manifold %s %s %d\n", argv[1], argv[2], resolution);

  obj.Process_Manifold(resolution);
  obj.SaveOBJ(argv[2]);
   return 0; 
}
