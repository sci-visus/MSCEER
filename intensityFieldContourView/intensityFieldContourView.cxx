/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkCamera.h"
#include "vtkConeSource.h"
#include "vtkContourFilter.h"
#include "vtkDataSet.h"
#include "vtkElevationFilter.h"
#include "vtkImageReader.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPiecewiseFunction.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkImageData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
//#include <vtkXMLPolyDataWriter.h>
#include "vtkQuadricDecimation.h"
#include <stdio.h>
//#include "vtkPLYWriter.h"
#include "vtkMarchingCubes.h"
#include "vtkFieldData.h"
#include "vtkProbeFilter.h"
#include "vtkPointData.h"
#include "vtkClipPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkProperty.h"
#include "vtkDoubleArray.h"
#include "vtkLookupTable.h"
#include "vtkPolyLine.h"
#include "vtkCellArray.h"
#include "vtkDecimatePro.h"
#include "vtkTriangle.h"
#include "vtkFeatureEdges.h"
#include "vtkKdTreePointLocator.h"
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>


#include <vtkGradientFilter.h>
#include <vtkImageRGBToHSV.h>
#include <vtkImageReader2.h>
#include <vtkImageReader2Factory.h>
#include <math.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageGradient.h>
#include <vtkImageShiftScale.h>
#include <vtkFieldDataToAttributeDataFilter.h>
#include <vtkChartXY.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkPlotPoints.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
#include "vtkImageGradientMagnitude.h"

#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkTable.h>

#include <vtkImageReslice.h>
#include <vtkTransform.h>
#include <vtkMarchingSquares.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageSlice.h>

#include <vtkKMeansStatistics.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkGraphicsFactory.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkMultiBlockDataSet.h>

#include <sstream>

#include <thread>

inline unsigned long LinearIndexFromCoordinate(int x, int y, int z, unsigned long X, unsigned long Y) {
  return x + y*X + z*(X*Y);//i*(N*M) + j*M + k; //index=y*WIDTH + x (2d)
}
int LinearIndexFromCoordinate(int x, int y, int z, int X, int Y) {
  return x + y*X + z*(X*Y);//i*(N*M) + j*M + k; //index=y*WIDTH + x (2d)
}
int print(int n){
    cout<<"place:" <<n <<std::endl;
    return 0;
}

void computeGradientField(vtkSmartPointer<vtkGradientFilter> gradientFilter, vtkSmartPointer<vtkImageData> subCubeGrad){
    // Compute the gradient of the Value

#if VTK_MAJOR_VERSION <= 5
    gradientFilter->SetInput(subCubeGrad);
#else
    gradientFilter->SetInputData(subCubeGrad);
#endif
    //SetDimensionality(3);
    //gradientFilter->ComputeRegularGridGradient()
    gradientFilter->Update(0);
    gradientFilter->Update(1);
}
void kMeansGradInensityField(vtkSmartPointer<vtkPoints> gradIntensityPoints, vtkSmartPointer<vtkKMeansStatistics> kMeansStatistics){
    // Get the points into the format needed for KMeans
    vtkSmartPointer<vtkTable> kMeansData =
        vtkSmartPointer<vtkTable>::New();

    for( int c = 0; c < 4; ++c )
      {
      std::stringstream colName;
      colName << "coord " << c;
      vtkSmartPointer<vtkFloatArray> gradintarray =
        vtkSmartPointer<vtkFloatArray>::New();
      gradintarray->SetNumberOfComponents(1);
      gradintarray->SetName( colName.str().c_str() );
      gradintarray->SetNumberOfTuples(gradIntensityPoints->GetNumberOfPoints());

      for( int r = 0; r < gradIntensityPoints->GetNumberOfPoints(); ++ r )
        {
        double p[3];
        gradIntensityPoints->GetPoint(r, p);

        gradintarray->SetValue( r, p[c] );
        }

      kMeansData->AddColumn( gradintarray );
      }

#if VTK_MAJOR_VERSION <= 5
    kMeansStatistics->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, kMeansData );
#else
    kMeansStatistics->SetInputData( vtkStatisticsAlgorithm::INPUT_DATA, kMeansData );
#endif
    kMeansStatistics->SetColumnStatus( kMeansData->GetColumnName( 0 ) , 1 );
    kMeansStatistics->SetColumnStatus( kMeansData->GetColumnName( 1 ) , 1 );
    kMeansStatistics->SetColumnStatus( kMeansData->GetColumnName( 2 ) , 1 );
    kMeansStatistics->SetColumnStatus( kMeansData->GetColumnName( 3 ) , 1 );
    //kMeansStatistics->SetColumnStatus( "Testing", 1 );
    kMeansStatistics->RequestSelectedColumns();
    kMeansStatistics->SetAssessOption( true );
    kMeansStatistics->SetDefaultNumberOfClusters( 4 );
    kMeansStatistics->Update() ;


}

void computeIntensityFieldIsosurface(vtkSmartPointer<vtkMarchingCubes> tIsoSurfProgress0, vtkSmartPointer<vtkImageData> subCubeIso, float isovalue){
    tIsoSurfProgress0->SetInputData(subCubeIso);
    tIsoSurfProgress0->SetValue(0, isovalue);
    tIsoSurfProgress0->Update(0);
}

// Using transormation matrix vmatrix compute slice vtkSlice of volume vtkdata
vtkImageData* compute_slice(vtkImageData *vtkdata, vtkMatrix4x4 *vmatrix, vtkAbstractTransform *vtktransform, vtkImageData *vtkSlice, vtkImageData *displaySlice) {
  //#ifndef USE_VTK
  //  print Visit not available!;
  //#else
  vtkSmartPointer<vtkImageReslice> vtkSlicer = vtkSmartPointer<vtkImageReslice>::New();
  vtkSlicer->SetOutputDimensionality(2);
  vtkSlicer->SetTransformInputSampling(1);
  // having auto crop on ensures nothing cropped however adds points of positive value to edge
  //vtkSlicer->SetAutoCropOutput(true);
  vtkSlicer->SetBorder(false);
  vtkSlicer->SetWrap(1);
  vtkSlicer->SetInterpolationModeToCubic();

  vtkSlicer->SetResliceAxes(vmatrix);

  vtkSlicer->SetInputData(vtkdata);
  vtkSlicer->Update();

  //vtkImageData *slice = vtkImageData::New();
  vtkSlice->DeepCopy(vtkSlicer->GetOutput());
  displaySlice->DeepCopy(vtkSlicer->GetOutput());
  //vtkSmartPointer<vtkTransform> vtkTransform = vtkSmartPointer<vtkTransform>::New();
  vtktransform = vtkSlicer -> GetResliceTransform();
  //vtkTransform = vtkSlicer -> ResliceTransform();
  return vtkSlice;

}

// return unit vector of v
std::vector<double> Unit(std::vector<double> v){
  float mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] *= 1.0/mag;
  v[1] *= 1.0/mag;
  v[2] *= 1.0/mag;
  return v;
}

// main function, takes command args:
// progressvariableraw heatreleaseraw x y z progressisoval heatisoval numberofcrosssections [reducesurfacetotarget#]
int main(int argc, char** argv) {


	// READ IN COMMAND LINE ARGUMENTS
    if (argc < 6) {
        printf("Command line args: File1 X Y Z (volume_fraction scale factor) (isovalue) (interactive_view) \n");
        exit(1);
	}

    float isovalue;
    sscanf(argv[6], "%f", &isovalue);
    int scaleFactor;
    int run_kMeans= 1;
    if(isovalue == 0.0){
        std::cout << "\n Using Kmeans clustering gradient field vs. intensity field centroids for isovalue" << std::endl;
    }else{
        run_kMeans=0;
    }
    sscanf(argv[5], "%d", &scaleFactor);

    std::cout << "\n Using scale factor: " << scaleFactor <<std::endl;

    int interactive;
    if(!argv[7])
        interactive = 0;
    else
        sscanf(argv[7], "%d", &interactive);



	int X,Y,Z;
    sscanf(argv[2], "%d", &X);
    sscanf(argv[3], "%d", &Y);
    sscanf(argv[4], "%d", &Z);
    //float clipvalue0;
    //sscanf(argv[7], "%f", &clipvalue0);


	// CONVERT the "c" value of heat release to the corresponding H2 value
    //isovalue = (isovalue * -0.0116013) + 0.011607718;

	// READ IN DATA
    std::cout << "\nReading In File: " << argv[1] << std::endl;
    vtkSmartPointer<vtkImageReader> tReaderProgress0 = vtkSmartPointer<vtkImageReader>::New();
	tReaderProgress0->SetFileName(argv[1]);
	tReaderProgress0->SetFileDimensionality(3);
	tReaderProgress0->SetDataExtent(0, X-1, 0, Y-1, 0, Z-1);
    tReaderProgress0->SetDataByteOrderToLittleEndian();
    tReaderProgress0->SetDataScalarTypeToFloat();
    tReaderProgress0->Update(0);

    vtkSmartPointer<vtkImageReader> tReaderProgress1 = vtkSmartPointer<vtkImageReader>::New();
    tReaderProgress1->SetFileName(argv[1]);
    tReaderProgress1->SetFileDimensionality(3);
    tReaderProgress1->SetDataExtent(0, X-1, 0, Y-1, 0, Z-1);
    tReaderProgress1->SetDataByteOrderToLittleEndian();
    tReaderProgress1->SetDataScalarTypeToFloat();
    tReaderProgress1->Update(0);

    vtkSmartPointer<vtkImageData> originalImage
            = vtkSmartPointer<vtkImageData>::New();
    originalImage->SetDimensions(X,Y,Z);
    originalImage->SetOrigin(0, 0, 0);
    originalImage->SetSpacing(1.0,1.0,1.0);
    originalImage = tReaderProgress0->GetOutput();

    vtkSmartPointer<vtkImageData> subCubeGrad
            = vtkSmartPointer<vtkImageData>::New();
    if(run_kMeans){
        subCubeGrad->SetOrigin(0.0, 0.0, 0.0);
        subCubeGrad->SetSpacing(1.0,1.0,1.0);
        subCubeGrad->SetDimensions((int)X/scaleFactor, (int)Y/scaleFactor, (int)Z/scaleFactor); //set dim to max-min of bbox
        subCubeGrad->AllocateScalars(VTK_FLOAT, 1);
    }

    vtkSmartPointer<vtkImageData> subCubeIso
            = vtkSmartPointer<vtkImageData>::New();
    subCubeIso->SetOrigin(0.0, 0.0, 0.0);
    subCubeIso->SetSpacing(1.0,1.0,1.0);
    subCubeIso->SetDimensions((int)X/scaleFactor, (int)Y/scaleFactor, (int)Z/scaleFactor); //set dim to max-min of bbox
    subCubeIso->AllocateScalars(VTK_FLOAT, 1);

    int subX = (int)X/scaleFactor;
    int subY = (int)Y/scaleFactor;
    int subZ = (int)Z/scaleFactor;
    int subRange = subX*subY*subZ;

    float* intensity_xyz
            = new float[subRange];

    int cropStartX = subX;
    int cropStartY = subY;
    int cropStartZ = subZ;
    if(scaleFactor>5 || scaleFactor<3){
        cropStartX = (int)X/3;
        cropStartY = (int)Y/3;
        cropStartZ = (int)Z/3;
    }

    for (int x = 0; x < (int)X/scaleFactor; x++){
      for (int y = 0; y < (int)Y/scaleFactor; y++){
        for (int z = 0; z < (int)Z/scaleFactor; z++){

          float* pixel1
                  = static_cast<float*>(tReaderProgress1->GetOutput()->GetScalarPointer(x+cropStartX,y+cropStartY,z+cropStartZ));
          if(run_kMeans==1){
            float* pixel0
                    = static_cast<float*>(tReaderProgress0->GetOutput()->GetScalarPointer(x+cropStartX,y+cropStartY,z+cropStartZ));

            float* pg
                    = static_cast<float*>(subCubeGrad->GetScalarPointer(x,y,z));
            *pg = *pixel0;
          }
          float* pi
                  = static_cast<float*>(subCubeIso->GetScalarPointer(x,y,z));


          *pi = *pixel1;
          //std::cout << pi[0] <<std::endl;
          //std::cout << pg[0] <<std::endl;
          intensity_xyz[LinearIndexFromCoordinate(x,y,z,subX, subY)]
                  = pi[0];
        }
      }
    }

    // Compute the gradient of the Value
    vtkSmartPointer<vtkGradientFilter> gradientFilter
            = vtkSmartPointer<vtkGradientFilter>::New();
    //dummy func to compile w/ no kmeans
    auto dummy = [](int hack){return hack+0;};
    std::thread gradientFieldThread;//(dummy, 0);
    if(run_kMeans==1){
        gradientFieldThread = std::thread(computeGradientField, gradientFilter, subCubeGrad);
    }
    // prepare slice matrix and plane information
    vtkImageData* vtkSlice = vtkImageData::New();
    vtkMatrix4x4* SliceMatrix;
    vtkAbstractTransform *vtkTransform;

    SliceMatrix = vtkMatrix4x4::New();
    SliceMatrix->Identity();

    std::vector<double> nz = {(double)X/(2.0), (double)Y/(2.0), (double)Z/2.0 };

    double oz = nz[2];

    //#x axis column
    SliceMatrix->SetElement(0, 0, 1);
    SliceMatrix->SetElement(1, 0, 0);
    SliceMatrix->SetElement(2, 0, 0);
    SliceMatrix->SetElement(3, 0, 0.0);

    //#y axis column
    SliceMatrix->SetElement(0, 1, 0);
    SliceMatrix->SetElement(1, 1, 1);
    SliceMatrix->SetElement(2, 1, 0);
    SliceMatrix->SetElement(3, 1, 0.0);

    //obtain normal of plane intersecting middle of
    // subcube in z component directed towards top
    //nz = Unit(nz);

    //#z axis column
    SliceMatrix->SetElement(0, 2, 0);//nz[0]);
    SliceMatrix->SetElement(1, 2, 0);//nz[1]);
    SliceMatrix->SetElement(2, 2, 1);//nz[2]);
    SliceMatrix->SetElement(3, 2, 0.0);

    SliceMatrix->SetElement(0, 3, nz[0]);
    SliceMatrix->SetElement(1, 3, nz[1]);
    SliceMatrix->SetElement(2, 3, nz[2]);
    SliceMatrix->SetElement(3, 3, 1.0);

    vtkImageData* sliceImage = vtkImageData::New();
    // Compute slice of cube given coordinated defining perpendicularly intersecting
    // plane to interior line
    //compute_slice(originalImage, SliceMatrix, vtkTransform, vtkSlice, sliceImage);

    std::thread sliceThread(compute_slice,originalImage, SliceMatrix, vtkTransform, vtkSlice, sliceImage);

    double xRange[2];
    if(run_kMeans){
        // merge gradient thread
        gradientFieldThread.join();

        // Extract the x component of the gradient
        vtkSmartPointer<vtkImageExtractComponents> extractXGradFilter
                = vtkSmartPointer<vtkImageExtractComponents>::New();
        extractXGradFilter->SetComponents(0);
        extractXGradFilter->SetInputConnection(gradientFilter->GetOutputPort());

    //    //double xCompRange[2];

        extractXGradFilter->Update(0);

        extractXGradFilter->GetOutput()->GetPointData()->GetScalars()->GetRange(xRange);


        // Extract the y component of the gradient
        vtkSmartPointer<vtkImageExtractComponents> extractYGradFilter =
          vtkSmartPointer<vtkImageExtractComponents>::New();
        extractYGradFilter->SetComponents(0);//1);
        extractYGradFilter->SetInputConnection(gradientFilter->GetOutputPort());

        //double yComp[2];
        extractYGradFilter->Update(0);
        vtkDataArray *yGrad = extractYGradFilter->GetOutput()->GetPointData()->GetScalars();
        //extractYGradFilter->GetOutput()->GetPointData()->GetScalars()->GetRange(yComp);

        vtkSmartPointer<vtkPoints> gradIntensityPoints =
            vtkSmartPointer<vtkPoints>::New();

        float* grad_xyz
                = new float[subRange];

        for (int x = 0; x < (int)X/scaleFactor; x++){
          for (int y = 0; y < (int)Y/scaleFactor; y++){
            for (int z = 0; z < (int)Z/scaleFactor; z++){

              float* pixel = static_cast<float*>(tReaderProgress0->GetOutput()->GetScalarPointer(x+cropStartX,y+cropStartY,z+cropStartZ));
              float* gradp = static_cast<float*>(gradientFilter->GetImageDataOutput()->GetScalarPointer(x,y,z));
              float* gradx = static_cast<float*>(extractXGradFilter->GetOutput()->GetScalarPointer(x,y,z));
              float* grady = static_cast<float*>(extractYGradFilter->GetOutput()->GetScalarPointer(x,y,z));

              int idx = LinearIndexFromCoordinate(x, y, z, subX , subY);

              float gx = gradx[0];
              if(gradx[0] != 0.0){
                  gx = log(gradx[0]);
              }
              float gy = grady[0];
              if(grady[0] != 0.0){
                  gy = log(grady[0]);
              }

              grad_xyz[idx] = gradx[0];

              gradIntensityPoints->InsertNextPoint(intensity_xyz[idx], gx, gy);

            }
          }
        }

        vtkSmartPointer<vtkKMeansStatistics> kMeansStatistics =
          vtkSmartPointer<vtkKMeansStatistics>::New();

        std::thread kmeansThread(kMeansGradInensityField, gradIntensityPoints, kMeansStatistics);

        //join kmeans thread
        kmeansThread.join();

        vtkMultiBlockDataSet* outputMetaDS = vtkMultiBlockDataSet::SafeDownCast( kMeansStatistics->GetOutputDataObject( vtkStatisticsAlgorithm::OUTPUT_MODEL ) );
        vtkSmartPointer<vtkTable> outputMeta = vtkTable::SafeDownCast( outputMetaDS->GetBlock( 0 ) );
        //vtkSmartPointer<vtkTable> outputMeta = vtkTable::SafeDownCast( outputMetaDS->GetBlock( 1 ) );
        vtkDoubleArray* coord0 = vtkDoubleArray::SafeDownCast(outputMeta->GetColumnByName("coord 0"));
        vtkDoubleArray* coord1 = vtkDoubleArray::SafeDownCast(outputMeta->GetColumnByName("coord 1"));
        vtkDoubleArray* coord2 = vtkDoubleArray::SafeDownCast(outputMeta->GetColumnByName("coord 2"));
        vtkDoubleArray* coord3 = vtkDoubleArray::SafeDownCast(outputMeta->GetColumnByName("coord 3"));

        for(unsigned int i = 0; i < coord0->GetNumberOfTuples(); ++i){
          std::cout <<"\nKMeans Cluster Centers: " << coord0->GetValue(i) << " " << coord1->GetValue(i) << " " << coord2->GetValue(i) << " " << coord3->GetValue(i) << std::endl;
        }
        std::cout << "\nobserved isoValues ("  << coord0->GetValue(0) << "," << coord0->GetValue(1) << "," << coord0->GetValue(2) << "," << coord0->GetValue(3) << ")" << std::endl;

        // %%%% Now compute isosurfaces of volume and slice of volume using observed
        // %%%% kmeans cluster centers starting with second to largest
        isovalue = (float)coord0->GetValue(1);
        std::cout << "\nUsing computed kMeans isoValue: " << isovalue << " " << "." << std::endl;

    } //END KMEANS


    // EXTRACT FIRST MARCHING CUBES SURFACE (of progress variable)
    vtkSmartPointer<vtkMarchingCubes> tIsoSurfProgress0 = vtkSmartPointer<vtkMarchingCubes>::New();
//    tIsoSurfProgress0->SetInputData(subCubeIso);
//    tIsoSurfProgress0->SetValue(0, isovalue);
//    tIsoSurfProgress0->Update(0);

    std::thread isoVolumeThread(computeIntensityFieldIsosurface, tIsoSurfProgress0, subCubeIso, isovalue);

    // join slice thread
    sliceThread.join();

    // EXTRACT FIRST MARCHING CUBES SURFACE (of progress variable)
    vtkSmartPointer<vtkMarchingSquares> tIsoSurfPlane = vtkSmartPointer<vtkMarchingSquares>::New();
    tIsoSurfPlane->SetInputData(vtkSlice);
    tIsoSurfPlane->SetValue(0, isovalue);
    tIsoSurfPlane->Update(0);

	// SET UP RENDERER
	vtkRenderer *renderer = 
		vtkRenderer::New();
	renderer->SetBackground(.1, .2, .1);
	vtkRenderWindow *renderWindow = 
		vtkRenderWindow::New();
	renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> interactor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    //join isovolumethread
    isoVolumeThread.join();
            vtkPolyDataMapper *mapper3 =
        vtkPolyDataMapper::New();
	mapper3->SetInputConnection(tIsoSurfProgress0->GetOutputPort());
	vtkActor *actor3 = 
		vtkActor::New();
    actor3->SetMapper(mapper3);
    renderer->AddActor(actor3);

    //render for slice isosurface
    vtkPolyDataMapper *mapperSlice
            = vtkPolyDataMapper::New();
    mapperSlice->SetInputConnection(tIsoSurfPlane->GetOutputPort());
    mapperSlice->ScalarVisibilityOff();
    vtkActor *actorSlice = vtkActor::New();
    actorSlice->SetMapper(mapperSlice);
    actorSlice->GetProperty()->SetColor(.9,0,0);
    renderer->AddActor(actorSlice);

    // scale actual slice
    // Scale the output (0,255)
    if(!run_kMeans){
        sliceImage->GetPointData()->GetScalars()->GetRange(xRange);
                //->GetOutput()->GetPointData()->GetScalars()->GetRange(xRange);
    }
    vtkSmartPointer<vtkImageShiftScale> shiftScaleX =
      vtkSmartPointer<vtkImageShiftScale>::New();
    shiftScaleX->SetOutputScalarTypeToUnsignedChar();
    shiftScaleX->SetScale(255.0 / (xRange[1]-1));
    shiftScaleX->SetInputData(sliceImage);//Connection(sliceImage->GetOutputPort());

    //render actual image slice
    vtkSmartPointer<vtkImageResliceMapper> imageResliceMapper = vtkSmartPointer<vtkImageResliceMapper>::New();
  #if VTK_MAJOR_VERSION <= 5
    imageResliceMapper->SetInputConnection(shiftScaleX->GetProducerPort());
  #else
    imageResliceMapper->SetInputConnection(shiftScaleX->GetOutputPort());
  #endif
    vtkSmartPointer<vtkImageSlice> imageSliceMapper = vtkSmartPointer<vtkImageSlice>::New();
    imageSliceMapper->SetMapper(imageResliceMapper);


      renderer->AddViewProp(imageSliceMapper);

      double viewpoint[3];
//    double tmpx, tmpz;
	
//    tmpx = viewpoint[0];
//    tmpz = viewpoint[2];
//    viewpoint[0] = tmpz;
//    viewpoint[2] = tmpx;
    double posFocal[3];
    posFocal[0] = X / (2.0*scaleFactor);
    posFocal[1] = Y / (2.0*scaleFactor);
    posFocal[2] = Z / (2.0*scaleFactor);

    double distance = std::max(X,Y);
    distance = std::max(distance,(double)Z);

    //double lookat[3];
    viewpoint[0] = 1.1 * X;
    viewpoint[1] = 1.1 * Y;
    viewpoint[2] = 1.1 * Z;

    // set camera position
    renderer->ResetCamera();
    vtkCamera *camera = renderer->GetActiveCamera();
    //camera->OrthogonalizeViewUp();
    //camera->Elevation(0.01);
    double d = camera->GetDistance();
    std::cout << "distance " << d << std::endl;
    //SetDistance(0.0);//viewpoint[0]);
    camera->SetFocalPoint(posFocal);
    //camera->SetPosition(viewpoint);
    camera->Zoom(2.9);
    camera->Azimuth(-65.0);
    //camera->Dolly(1.3);

    renderWindow->Render();

    vtkSmartPointer<vtkWindowToImageFilter> wti = vtkSmartPointer<vtkWindowToImageFilter>::New();

    wti->SetInput(renderWindow); //Connection(renderWindow->GetOffScreenRendering());
    wti->SetInputBufferTypeToRGBA();

    wti->ReadFrontBufferOff();
    wti->Update(0);

    std::cout << "writing png file" << std::endl;

    auto png = vtkSmartPointer<vtkImageWriter>::New();
    png = vtkSmartPointer<vtkPNGWriter>::New();
    char tmpout[1024];
    sprintf(tmpout, "%s_isosurface_%f.png", argv[1],isovalue);
    printf("printing to file: %s \n", tmpout);
    png->SetFileName(tmpout);
    png->SetInputConnection(wti->GetOutputPort());
	png->Write();

	//return 0;
    if(interactive==1)
        interactor->Start();

	return 0;
}
