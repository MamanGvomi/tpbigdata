#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>

#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkShortArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkGraphicsFactory.h>
//#include <vtkImagingFactory.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkProperty.h>
#include <vtkVectorNorm.h>
#include <time.h>
#include <math.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>


#include "config.h"
#include "unistd.h"
#include <mpi.h>

#define BIG


/*
#define FICHIER MY_MESHES_PATH "/Frog.raw_256_256_44_CHAR.raw"

int gridSize = 256;
int YgridSize = 256;
int ZgridSize = 44;

#define CHAR

int startexploreval=13;
int endexploreval=20;

*/


//#define FICHIER MY_MESHES_PATH "/Mystere5_2048_2048_756_SHORT.raw"
//
// int gridSize = 2048;
// int YgridSize = 2048;
// int ZgridSize = 756;
//
// #define SHORT
//
//int startexploreval=54000;
//int endexploreval=55000;



//#define FICHIER MY_MESHES_PATH "/Mystere4_512_512_322_SHORT.raw"
//
// int gridSize = 512;
// int YgridSize = 512;
// int ZgridSize = 322;
//
// #define SHORT
//
//int startexploreval=50000;
//int endexploreval=60000;

//#define FICHIER MY_MESHES_PATH "/Mystere6_1118_2046_694_CHAR.raw"
//
//int gridSize = 1118;
//int YgridSize = 2046;
//int ZgridSize = 694;
//
//#define CHAR
//
//int startexploreval=15;
//int endexploreval=50;

#define FICHIER MY_MESHES_PATH "/Mystere8_2048_2048_2048_CHAR.raw"

int gridSize = 2048;
int YgridSize = 2048;
int ZgridSize = 2048;

#define CHAR

int startexploreval = 30;
int endexploreval = 80;


//#define FICHIER MY_MESHES_PATH "/Mystere1_512_512_134_SHORT.raw"
//int gridSize = 512;
//int YgridSize = 512;
//int ZgridSize = 134;
//
//#define SHORT
//
//int startexploreval = 30000;
//int endexploreval = 45000;

//#define FICHIER MY_MESHES_PATH "/Mystere2_512_400_512_SHORT.raw"
//int gridSize = 512;
//int YgridSize = 400;
//int ZgridSize = 512;
//
//#define SHORT
//
//int startexploreval = 21000;
//int endexploreval = 25000;

/*
#define FICHIER MY_MESHES_PATH "/Mystere4_512_512_322_SHORT.raw"
int gridSize = 512;
int YgridSize = 512;
int ZgridSize = 322;

 #define SHORT

 int startexploreval=1;
 int endexploreval=65000;
//*/


const char *location = FICHIER;


int winSize = 500;

int numPasses = 5;
int nbimages = 5;


const char *prefix = "";


int passNum = 0;

using std::cerr;
using std::endl;

// Function prototypes
vtkRectilinearGrid *ReadGrid(int zStart, int zEnd);

vtkRectilinearGrid *ParallelReadGrid(void);


void WriteImage(const char *name, const float *rgba, int width, int height);

bool CompositeImage(float *rgba_in, float *zbuffer, float *rgba_out,
                    int image_width, int image_height);

bool ComposeImageZbuffer(float *rgba_out, float *zbuffer, int image_width, int image_height);

int parRank = 0;
int parSize = 1;

int main(int argc, char *argv[]) {

    bool once = false;
    bool parallele = true;

    // MPI setup
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &parRank);
    MPI_Comm_size(MPI_COMM_WORLD, &parSize);

    srand(time(NULL));
//    int zStart = 0;
//    int zEnd = ZgridSize;

    int npixels = winSize * winSize;
    vtkRectilinearGrid *reader = NULL;
    if (parallele) {
        vtkRectilinearGrid *rg = ParallelReadGrid();

    }

    //LookupTable
    vtkLookupTable *lut = vtkLookupTable::New();
    lut->SetHueRange(0.1, 0.0);
    lut->SetSaturationRange(0.0, 1.0);
    lut->SetValueRange(1.0, 255.0);
    lut->SetNumberOfColors(100);
    lut->Build();

    // Renderer
    vtkRenderer *ren = vtkRenderer::New();
//    double bounds[6] = {0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001};
//    ren->ResetCamera(bounds);

    //Windows
    vtkRenderWindow *renwin = vtkRenderWindow::New();
    renwin->SetSize(winSize, winSize);
    renwin->AddRenderer(ren);

    // Filtre de contour
    vtkContourFilter *cf = vtkContourFilter::New();
    cf->SetNumberOfContours(1);
    // cf->SetValue(0, 20.0);
    int valcont = startexploreval;
    cf->SetValue(1, valcont);
    if (once) {
        //On lis tous le fichier
        int zStart = 0;
        int zEnd = ZgridSize;
        reader = (parallele) ? ParallelReadGrid() : ReadGrid(zStart, zEnd);
    }
    //Sinon y'a rien pour le moment
    cf->SetInputData(reader);


    // Ici je crois qu'on s'assure que l'image s'affiche dans la windows
    int maxsize = std::max(gridSize, std::max(YgridSize, ZgridSize));
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Scale(gridSize / (float) maxsize, YgridSize / (float) maxsize, ZgridSize / (float) maxsize);
//    transform->RotateY(+40);
//    transform->RotateX(-20);
    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputConnection(cf->GetOutputPort());
//    transformFilter->SetInputData(reader);
    transformFilter->SetTransform(transform);


    // Mapper
    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
//    mapper->SetInputData(reader);
    mapper->SetInputConnection(transformFilter->GetOutputPort());

    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);

    float startval = valcont;
    float endval = endexploreval;


    mapper->SetScalarRange(startexploreval, endexploreval);
    mapper->SetLookupTable(lut);

    ren->AddActor(actor);
    ren->SetViewport(0, 0, 1, 1);

    vtkCamera *cam;

    // vtkCamera *cam = ren->GetActiveCamera();
//            if(once){
//                once=false;
//                cam= ren->GetActiveCamera();
//             cam->SetFocalPoint(0.5, 0.5, 0.5);
//                cam->SetPosition(0., .0, 3.);
//                cam->SetViewUp(0., -1.0, 0.0);
//                bounds=ren->ComputeVisiblePropBounds();
//                std::cout<<"bounds : "<<bounds[0]<<" ; "<<bounds[1]<<" ; "<<bounds[2]<<" ; "<<bounds[3]<<" ; "<<bounds[4]<<" ; "<<bounds[5]<<" ; "<<endl<<fflush;
//            }
//            else{ren->SetActiveCamera(cam);
//
//            }

    cam = ren->GetActiveCamera();
    cam->SetFocalPoint(0.5, 0.5, 0.5);
    cam->SetPosition(-1., 0.0, 3.);
    cam->SetViewUp(0, 1.0, 0.0);
//    cam->Roll(-180);
//    ren->ResetCamera(bounds);

//    cam->Pitch();
//    bounds=ren->ComputeVisiblePropBounds();
    // std::cout<<"bounds : "<<bounds[0]<<" ; "<<bounds[1]<<" ; "<<bounds[2]<<" ; "<<bounds[3]<<" ; "<<bounds[4]<<" ; "<<bounds[5]<<" ; "<<endl<<fflush;

    //       ren->SetActiveCamera (cam);
    //           cam= ren->GetActiveCamera();
    //   cam->SetPosition(0.5, 3.0, 0.5);
    //  cam->SetViewUp(0., 0.0, 1.0);
    //      cam->SetFocalPoint(0., 0.0, 0.);

//

    //*

    if (once) {
        vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
        iren->SetRenderWindow(renwin);
        renwin->Render();
        iren->Start();
    } else {
        //Les buffers locaux
//        renwin->Render();
        float *rgba = new float[4 * winSize * winSize];
        float *zbuffer = new float[winSize * winSize];
        float *auxrgba = new float[4 * winSize * winSize]; // rgba pour final
        float *auxzbuffer = new float[4 * winSize * winSize]; // zbuffer pour final
        bool transparence = false;


        //Init des buffers pour l'image final
        for (int i = 0; i < winSize * winSize; i++) {
            auxzbuffer[i] = 1.0;
            auxrgba[i * 4] = 0;
            auxrgba[i * 4 + 1] = 0;
            auxrgba[i * 4 + 2] = 1;
            auxrgba[i * 4 + 3] = 0;
        }
        //Je ne gère pas les restes.... :'(
        if (!parallele) {

            int step = (ZgridSize / nbimages);
            for (numPasses = 0; numPasses < nbimages; numPasses++) {

                //La zone de travail
                int zEnd = (numPasses + 1) * step;
                int zStart = numPasses * step;

                //Sinon ca bug a la fin (je sais pas trop pourquoi honnetement)
                if (numPasses == nbimages - 1) {
                    zEnd = (numPasses + 1) * step - 1;
                }

                //On remplis
                reader = (parallele) ? ParallelReadGrid() : ReadGrid(zStart, zEnd);
                //Le filtre conteur le lis
                cf->SetInputData(reader);
                //On le delete (il a fait son taf)
                reader->Delete();

                // Force an update and set the parallel rank as the active scalars.
                cf->Update();
//            double * bounds=ren->ComputeVisiblePropBounds();
//            ren->ResetCamera(bounds);
                renwin->Render();

                rgba = renwin->GetRGBAPixelData(0, 0, winSize - 1, winSize - 1, 1);
                zbuffer = renwin->GetZbufferData(0, 0, winSize - 1, winSize - 1);
                char name[128];
                sprintf(name, "image%d.png", numPasses);
                WriteImage(name, rgba, winSize, winSize);


                float *new_rgba = new float[4 * winSize * winSize];
                bool didComposite = (parallele) ? CompositeImage(rgba, zbuffer, new_rgba, winSize, winSize)
                                                : ComposeImageZbuffer(new_rgba, zbuffer, winSize, winSize);

                char namez[128];
                sprintf(namez, "imageZ%d.png", numPasses);
                WriteImage(namez, new_rgba, winSize, winSize);
                // Pour l'image finale
                for (int i = 0; i < winSize * winSize; i++) {
                    if (!transparence) {
                        if (auxzbuffer[i] >= zbuffer[i]) {
                            auxzbuffer[i] = zbuffer[i];
                            auxrgba[i * 4] = rgba[i * 4];
                            auxrgba[i * 4 + 1] = rgba[i * 4 + 1];
                            auxrgba[i * 4 + 2] = rgba[i * 4 + 2];
                            auxrgba[i * 4 + 3] = rgba[i * 4 + 3];
                        }
                    } else {
                        auxrgba[i * 4] = (auxrgba[i * 4] + rgba[i * 4]) / 2 + 510;
                        auxrgba[i * 4 + 1] = (auxrgba[i * 4 + 1] + rgba[i * 4 + 1]) / 2 + 510;
                        auxrgba[i * 4 + 2] = (auxrgba[i * 4 + 2] + rgba[i * 4 + 2]) / 2 + 510;
                        auxrgba[i * 4 + 3] = (auxrgba[i * 4 + 3] + rgba[i * 4 + 3]) / 2 + 510;
                    }

                }

                free(rgba);
                free(zbuffer);
                free(new_rgba);
            }

        } else {
            reader = ParallelReadGrid();
            //Le filtre conteur le lis
            cf->SetInputData(reader);
            //On le delete (il a fait son taf)
            reader->Delete();
            // Force an update and set the parallel rank as the active scalars.
            cf->Update();
//            double * bounds=ren->ComputeVisiblePropBounds();
//            ren->ResetCamera(bounds);
            renwin->Render();
            rgba = renwin->GetRGBAPixelData(0, 0, (winSize - 1), (winSize - 1), 1);
            zbuffer = renwin->GetZbufferData(0, 0, (winSize - 1), (winSize - 1));


            float *new_rgba = new float[4 * winSize * winSize];
            bool didComposite = CompositeImage(rgba, zbuffer, new_rgba, winSize, winSize);
            if (didComposite) {
                if (parRank == 0) {
                    WriteImage("final_image_MPI.png", new_rgba, winSize, winSize);
                }

                char name[128];
                sprintf(name, "img_MPI%d.png", parRank);
                WriteImage(name, rgba, winSize, winSize);

            }
        }
//        WriteImage((transparence) ? "final_image_transparence.png" : "final_image.png", auxrgba, winSize, winSize);
    }
    //*/

    reader->Delete();
    mapper->Delete();
    cf->Delete();
    //          lut->Delete();
    ren->RemoveActor(actor);
    actor->Delete();

    ren->Delete();
    renwin->Delete();


}



// You should not need to modify these routines.

vtkRectilinearGrid *
ReadGrid(int zStart, int zEnd) {
    int i;

    /*   if (zStart < 0 || zEnd < 0 || zStart >= gridSize || zEnd >= gridSize || zStart > zEnd)
     {
     cerr << prefix << "Invalid range: " << zStart << "-" << zEnd << endl;
     return NULL;
     }
     */
    ifstream ifile(location);
    if (ifile.fail()) {
        cerr << prefix << "Unable to open file: " << location << "!!" << endl;
        throw std::runtime_error("can't find the file!! Check the name and the path of this file? ");
    }

    cerr << prefix << "Reading from " << zStart << " to " << zEnd << endl;

    vtkRectilinearGrid *rg = vtkRectilinearGrid::New();
    vtkFloatArray *X = vtkFloatArray::New();
    X->SetNumberOfTuples(gridSize);
    for (i = 0; i < gridSize; i++)
        X->SetTuple1(i, i / (gridSize - 1.0));
    rg->SetXCoordinates(X);
    X->Delete();
    vtkFloatArray *Y = vtkFloatArray::New();
    Y->SetNumberOfTuples(YgridSize);
    for (i = 0; i < YgridSize; i++)
        Y->SetTuple1(i, i / (YgridSize - 1.0));
    rg->SetYCoordinates(Y);
    Y->Delete();
    vtkFloatArray *Z = vtkFloatArray::New();
    int numSlicesToRead = zEnd - zStart + 1;
    Z->SetNumberOfTuples(numSlicesToRead);
    for (i = zStart; i <= zEnd; i++)
        Z->SetTuple1(i - zStart, i / (ZgridSize - 1.0));
    rg->SetZCoordinates(Z);
    Z->Delete();

    rg->SetDimensions(gridSize, YgridSize, numSlicesToRead);

    unsigned int valuesPerSlice = gridSize * YgridSize;

#if defined(SHORT)
    unsigned int bytesPerSlice = sizeof(short) * valuesPerSlice;

#elif defined(CHAR)
    unsigned int bytesPerSlice = sizeof(char) * valuesPerSlice;

#elif  defined(FLOAT)
    unsigned int bytesPerSlice   = sizeof(float)*valuesPerSlice;

#else
#error Unsupported choice setting
#endif


#if defined(SMALL)
    unsigned int offset          = (unsigned int)zStart * (unsigned int)bytesPerSlice;
   unsigned int bytesToRead     = bytesPerSlice*numSlicesToRead;
   unsigned int valuesToRead    = valuesPerSlice*numSlicesToRead;
#elif defined(BIG)
    unsigned long long offset = (unsigned long long) zStart * bytesPerSlice;
    unsigned long long bytesToRead = (unsigned long long) bytesPerSlice * numSlicesToRead;
    unsigned int valuesToRead = (unsigned int) valuesPerSlice * numSlicesToRead;
#else
#error Unsupported choice setting
#endif


#if defined(SHORT)
    vtkUnsignedShortArray *scalars = vtkUnsignedShortArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned short *arr = scalars->GetPointer(0);

#elif defined(CHAR)
    vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned char *arr = scalars->GetPointer(0);

#elif  defined(FLOAT)
    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    float *arr = scalars->GetPointer(0);
#else
#error Unsupported choice setting
#endif


    ifile.seekg(offset, ios::beg);
    ifile.read((char *) arr, bytesToRead);
    ifile.close();

    int min = +2147483647;
    int max = 0;


    for (int i = 0; i < valuesToRead; i++) {
        if (min > (scalars->GetPointer(0))[i]) min = (scalars->GetPointer(0))[i];
        if (max < (scalars->GetPointer(0))[i]) max = (scalars->GetPointer(0))[i];

        if (rand() % (valuesToRead / 20) == 0) {
#if defined(SHORT)
            std::cout << (unsigned short) (scalars->GetPointer(0))[i] << " ";

#elif defined(CHAR)
            std::cout << +(unsigned char) (scalars->GetPointer(0))[i] << " ";

#elif  defined(FLOAT)
            std::cout<<(float)(scalars->GetPointer(0))[i]<<" ";

#else
#error Unsupported choice setting
#endif


        }
    }


    std::cout << endl << fflush;
    std::cout << "min value read: " << min << endl << fflush;
    std::cout << "max value read: " << max << endl << fflush;


    scalars->SetName("entropy");
    rg->GetPointData()->AddArray(scalars);
    scalars->Delete();

    vtkFloatArray *pr = vtkFloatArray::New();
    pr->SetNumberOfTuples(valuesToRead);
    for (int i = 0; i < valuesToRead; i++)
        pr->SetTuple1(i, passNum);

    pr->SetName("pass_num");
    rg->GetPointData()->AddArray(pr);
    pr->Delete();

    rg->GetPointData()->SetActiveScalars("entropy");

    cerr << prefix << " Done reading" << endl;
    return rg;
}


void
WriteImage(const char *name, const float *rgba, int width, int height) {
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
#if VTK_MAJOR_VERSION <= 5
    img->SetNumberOfScalarComponents(3);
    img->SetScalarTypeToUnsignedChar();
#else
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
#endif

    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++) {
            unsigned char *ptr = (unsigned char *) img->GetScalarPointer(i, j, 0);
            int idx = j * width + i;
            ptr[0] = (unsigned char) (255 * rgba[4 * idx + 0]);
            ptr[1] = (unsigned char) (255 * rgba[4 * idx + 1]);
            ptr[2] = (unsigned char) (255 * rgba[4 * idx + 2]);
        }


    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(name);
    writer->Write();

    img->Delete();
    writer->Delete();
}


bool ComposeImageZbuffer(float *rgba_out, float *zbuffer, int image_width, int image_height) {
    int npixels = image_width * image_height;

    float min = 1;
    float max = 0;
    for (int i = 0; i < npixels; i++) {
        if (zbuffer[i] < min) min = zbuffer[i];
        if (zbuffer[i] > max) max = zbuffer[i];

    }
    std::cout << "min:" << min;
    std::cout << "max:" << max << "  ";

    float coef = 1.0 / ((max - min));

    std::cout << "coef:" << coef << "  ";

    for (int i = 0; i < npixels; i++) {

        rgba_out[i * 4] = (zbuffer[i] == 1.0 ? 0 : 1 - coef * (zbuffer[i] - min));
        rgba_out[i * 4 + 1] = (zbuffer[i] == 1.0 ? 0 : 1 - coef * (zbuffer[i] - min));
        rgba_out[i * 4 + 2] = (zbuffer[i] == 1.0 ? 0 : 1 - coef * (zbuffer[i] - min));
        rgba_out[i * 4 + 3] = 0.0;
    }


    return false;
}


vtkRectilinearGrid *ParallelReadGrid(void) {
//    int zStart = (gridSize / 1.7) - 1;
//    int zEnd = (gridSize / 1.5) - 1;
    int zStart = (ZgridSize / parSize) * parRank;
    int zEnd = ((ZgridSize / parSize)) * (parRank + 1);
    if (parRank == parSize - 1) {
        zEnd = ((ZgridSize / parSize)) * (parRank + 1) - 1;
    }
//    int zStart = 0;
//    int zEnd = (gridSize)-1;

    // if you don't want any data for this processor, just do this:
    //    return NULL;

    return ReadGrid(zStart, zEnd);
}

bool CompositeImage(float *rgba_in, float *zbuffer, float *rgba_out,
                    int image_width, int image_height) {
    int npixels = image_width * image_height;
    float *zbufferMin = new float[npixels];
    float *rgba_tmp = new float[4 * npixels];

    MPI_Allreduce(zbuffer, zbufferMin, npixels, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    for (int i = 0; i < npixels; i++) {
        // Sans transparance
        if (zbufferMin[i] == zbuffer[i]) {
            rgba_tmp[i * 4] = rgba_in[i * 4];
            rgba_tmp[i * 4 + 1] = rgba_in[i * 4 + 1];
            rgba_tmp[i * 4 + 2] = rgba_in[i * 4 + 2];
            rgba_tmp[i * 4 + 3] = rgba_in[i * 4 + 3];
        }
        // Avec transparence
//            rgba_tmp[i * 4] = (rgba_tmp[i * 4] + rgba_in[i * 4]) / 2;
//            rgba_tmp[i * 4 + 1] =  (rgba_tmp[i * 4 + 1] + rgba_in[i * 4 + 1]) / 2;
//            rgba_tmp[i * 4 + 2] =  (rgba_tmp[i * 4 + 2] + rgba_in[i * 4 + 2]) / 2;
//            rgba_tmp[i * 4 + 3] =  (rgba_tmp[i * 4 + 3] + rgba_in[i * 4 + 3]) / 2;

    }
    MPI_Reduce(rgba_tmp, rgba_out, 4 * npixels, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    return true;
}