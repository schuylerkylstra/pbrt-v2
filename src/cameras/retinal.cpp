#include "stdafx.h"
#include "cameras/retinal.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include <iostream>

// RetinalCamera Definitions
RetinalCamera::RetinalCamera(const AnimatedTransform &cam2world, const float screenWindow[4],
					float sopen, float sclose, Film *film, float pupilradius, Point nodalP,
					const float *zVals, const Normal *Normals , float rIndex , float lensRadius) 

	: ProjectiveCamera(cam2world, Orthographic(0., 1.), screenWindow, sopen, sclose, 0.0, 1e30f, film)
{
	// compute differential changes in origin for retinal camera rays
	dxCamera = RasterToCamera(Vector(1, 0, 0));
	dyCamera = RasterToCamera(Vector(0, 1, 0));

    numCones    = 131072;
    xResolution = 512;
    yResolution = 256;
    NodalPoint  = nodalP;
    PupilRadius = pupilradius;
    Cornea      = Lens(rIndex, lensRadius, NodalPoint);
    zFlag       = false;
    nFlag       = false;


    // set the z values and normals of the retina
    if (zVals) 
    {
        zFlag = true;
        zValues = new float[numCones];
        memcpy(zValues, zVals, numCones*sizeof(float));
    }
    if (Normals) 
    {
        nFlag = true;
        rNormals = new Normal[numCones];
        memcpy(rNormals, Normals, numCones*sizeof(Normal));
    }
}

RetinalCamera::~RetinalCamera()
{
    if (zFlag)
    {
        delete[] zValues;
    }
    if (nFlag)
    {
        delete[] rNormals;
    }
}


float RetinalCamera::GenerateRay( const CameraSample &sample, Ray * ray) const
{
	// The new Origin =  (sample.imageX, sample.imageY, zVals[sample.imageX + sample.imageY * rowL])
	// the time should be set to 0
    // std::cout << "RAY\n";
    int row = int( sample.imageY );
    int col = int( sample.imageX );
    float z = zValues[ row * xResolution + col ];

    Point Pras(sample.imageX, sample.imageY, z);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);

    Vector Dras(NodalPoint - Pras);
    Vector Dcam;
    RasterToCamera(Dras, &Dcam);

    *ray = Ray(Pcamera, Normalize(Dcam), 0., INFINITY);

    if(PupilRadius > 0.f)
    {
        // compute offset
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= PupilRadius;
        lensV *= PupilRadius;

        // compute new ray
        Point Offset( 256.0 + lensU, 128.0 + lensV, 0.0);
        Vector ODir( Offset - Pras );
        *ray = Ray(Pras, ODir, 0., INFINITY);

        Cornea.ApplySnellLaw( ray );

        RasterToCamera( ray->o, &(ray->o) );
        RasterToCamera( ray->d, &(ray->d) );

        ray->d = Normalize( ray->d );

    }
    
	// This produces a vector from the retina to the nodal point
	// All points that pass through the nodal point are unaltered
	// need to alter this to include pupilradius information
	ray->time = sample.time;
	CameraToWorld(*ray, ray);
	return 1.f;
}


float RetinalCamera::GenerateRayDifferential(const CameraSample &sample,
        RayDifferential *ray) const {

    // Generate raster and camera samples
    // std::cout << "RAYDIFF\n";
    int row = int( sample.imageY );
    int col = int( sample.imageX );
    float z = zValues[ row * xResolution + col ];

    // convert from Raster to camera space
    Point Pras(sample.imageX, sample.imageY, z);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);

    Vector Dras(NodalPoint - Pras);
    Vector Dcam;
    RasterToCamera(Dras, &Dcam);

    // Generate ray differential
    *ray = RayDifferential( Pcamera, Normalize(Dcam), 0., INFINITY );

    ray->time = sample.time;

    CameraSample sshift = sample;


    Point rxNew(++sshift.imageX, sshift.imageY, zValues[ row * xResolution + col + 1 ] );
    Point ryNew(--sshift.imageX, ++sshift.imageY, zValues[ (row + 1) * xResolution + col ] );

    Vector rxDir, ryDir;

    if(PupilRadius > 0)
    {
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= PupilRadius;
        lensV *= PupilRadius;

        Point NodalOffset(256.0 + lensU, 128.0 + lensV, 0.0);
        Vector Odir( NodalOffset - Pras );

        rxDir = NodalPoint - rxNew;
        ryDir = NodalPoint - ryNew;

        Ray xRay(rxNew, rxDir, 0., INFINITY);
        Ray yRay(ryNew, ryDir, 0., INFINITY);
        Ray r(Pras, Odir, 0., INFINITY);

        Cornea.ApplySnellLaw( &r );
        Cornea.ApplySnellLaw( &xRay );
        Cornea.ApplySnellLaw( &yRay );

        // std::cout << r.o.x << ", " << r.o.y << ", " << r.o.z << "\n";
        // std::cout << r.d.x << ", " << r.d.y << ", " << r.d.z << "\n";
        Vector ODcam;
        Vector RDcam;
        Point  Impact;
        RasterToCamera( Odir, &ODcam );
        RasterToCamera( r.d, &RDcam );
        RasterToCamera( r.o - r.d*5, &Impact );
        // std::cout << Impact.x << ", " << Impact.y << ", " << Impact.z << "\n";
        // std::cout << RDcam.x << ", " << RDcam.y << ", " << RDcam.z << "\n";
        // *ray = RayDifferential( Pcamera, Normalize(ODcam), 0., INFINITY );
        *ray = RayDifferential( Impact, Normalize(RDcam), 0., INFINITY );

        rxNew = xRay.o;
        ryNew = yRay.o;

        rxDir = xRay.d;
        ryDir = yRay.d;
    }
    else
    {
        rxDir = NodalPoint - rxNew;
        ryDir = NodalPoint - ryNew;
    }
    
    

    RasterToCamera(rxNew, &(ray->rxOrigin) );
    RasterToCamera(ryNew, &(ray->ryOrigin) );
    RasterToCamera(rxDir, &(ray->rxDirection) );
    RasterToCamera(ryDir, &(ray->ryDirection) );

    ray->hasDifferentials = true;
    CameraToWorld(*ray, ray);
    return 1.f;
}


RetinalCamera* CreateRetinalCamera(const ParamSet &params, const AnimatedTransform &cam2world, Film *film)
{
	// Extract common camera parameters from _ParamSet_
    float shutteropen  = params.FindOneFloat("shutteropen", 0.f);
    float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        swap(shutterclose, shutteropen);
    }

    float frame           = params.FindOneFloat("frameaspectratio", float(film->xResolution)/float(film->yResolution));
    float lensradius      = params.FindOneFloat("lensradius", 5.5);
    Point nodalpoint      = params.FindOnePoint("nodalpoint", Point(256.0, 128.0, 0.0) );
    float pupilradius     = params.FindOneFloat("pupilradius", 0.f);
    float refractiveIndex = params.FindOneFloat("refractiveindex", 1.33);
    float screen[4];

    if (frame > 1.f) {
        screen[0] = -frame;
        screen[1] =  frame;
        screen[2] = -1.f;
        screen[3] =  1.f;
    }

    else {
        screen[0] = -1.f;
        screen[1] =  1.f;
        screen[2] = -1.f / frame;
        screen[3] =  1.f / frame;
    }

    int swi, zvi, rni;

    const float *sw         = params.FindFloat("screenwindow", &swi);
    const float *zVals      = params.FindFloat("zvalues", &zvi);
    const Normal *Normals   = params.FindNormal("normals", &rni);


    if (sw && swi == 4)
    {
        memcpy( screen, sw, 4*sizeof(float) );
    }
    
    return new RetinalCamera(cam2world, screen, shutteropen, shutterclose,
        film, pupilradius, nodalpoint, zVals, Normals, refractiveIndex, lensradius);
}
