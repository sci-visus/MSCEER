
#include <vector>
#include <set>

#include <vtkActor.h>
#include <vtkMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYWriter.h>

#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkInteractionStyle)

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkVoxelModeller.h>
#include <vtkSphereSource.h>
#include <vtkDICOMImageReader.h>

vtkSmartPointer<vtkImageData> getImageData(char* data, int* dims){
  vtkSmartPointer<vtkImageData> volume =
  vtkSmartPointer<vtkImageData>::New();
  
  // Specify the size of the image data
  volume->SetDimensions(dims[0], dims[1], dims[2]);
#if VTK_MAJOR_VERSION <= 5
  volume->SetNumberOfScalarComponents(1);
  volume->SetScalarTypeToFloat();
#else
  volume->AllocateScalars(VTK_DOUBLE,1);
#endif
  
//  std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
//  
//  std::cout << "Number of points: " << volume->GetNumberOfPoints() << std::endl;
//  std::cout << "Number of cells: " << volume->GetNumberOfCells() << std::endl;
  for (int z = 0; z < dims[2]; z++)
  {
    for (int y = 0; y < dims[1]; y++)
    {
      for (int x = 0; x < dims[0]; x++)
      {
        uint32_t idx = x + y * dims[0] + z*dims[0]*dims[1];
        double* pixel = static_cast<double*>(volume->GetScalarPointer(x,y,z));
        float pv = ((float*)data)[idx];
        pixel[0] = ((float*)data)[idx]; // TODO use template for the type
      }
    }
  }
  
  return volume;
  
}

int main(int argc, char** argv) {

	int X, Y, Z;
	std::string filename;
	if (argc < 5) { printf("Usage: X Y Z filename [isovalue]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	float isovalue = 0.0;
  if(argc > 5)
    sscanf(argv[5], "%f", &isovalue);
  
  int dims[3]={X,Y,Z};
  size_t mysize=dims[0]*dims[1]*dims[2]*sizeof(float);
  
  char* data = new char[mysize];
  
  std::ifstream file (filename, ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    std::streampos size;
    size = file.tellg();
    file.seekg (0, ios::beg);
    file.read (data, size);
    file.close();
    
    printf("read file done\n");
  }
  
  vtkSmartPointer<vtkImageData> volume = getImageData(data, dims);
  
  vtkSmartPointer<vtkMarchingCubes> surface =
  vtkSmartPointer<vtkMarchingCubes>::New();
  
#if VTK_MAJOR_VERSION <= 5
  surface->SetInput(volume);
#else
  surface->SetInputData(volume);
#endif
  surface->ComputeNormalsOn();
  surface->SetValue(0, isovalue);
  
  printf("computing isosurface for isovalue %f...\n", isovalue);

//  surface->Update();
  
//  vtkSmartPointer<vtkRenderer> renderer =
//  vtkSmartPointer<vtkRenderer>::New();
//  renderer->SetBackground(.1, .2, .3);
//  
//  vtkSmartPointer<vtkRenderWindow> renderWindow =
//  vtkSmartPointer<vtkRenderWindow>::New();
//  renderWindow->AddRenderer(renderer);
//  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
//  vtkSmartPointer<vtkRenderWindowInteractor>::New();
//  interactor->SetRenderWindow(renderWindow);
//  
//  vtkSmartPointer<vtkPolyDataMapper> mapper =
//  vtkSmartPointer<vtkPolyDataMapper>::New();
//  mapper->SetInputConnection(surface->GetOutputPort());
//  mapper->ScalarVisibilityOff();
//  
//  vtkSmartPointer<vtkActor> actor =
//  vtkSmartPointer<vtkActor>::New();
//  actor->SetMapper(mapper);
//  
//  renderer->AddActor(actor);
//  
//  renderWindow->Render();
//  interactor->Start();

  vtkSmartPointer<vtkPLYWriter> plyWriter =
  vtkSmartPointer<vtkPLYWriter>::New();
  plyWriter->SetInputConnection(surface->GetOutputPort());
  plyWriter->SetFileName((filename+".isosurface.ply").c_str());
  plyWriter->Write();
  
//  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
//  writer->SetInputData(marched);
//  writer->SetFileName("surface.vtk");
//  writer->SetFileTypeToASCII();
//  writer->Write();

	return 0;

}


