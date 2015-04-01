#include "stdafx.h"
#include "cameras/retinal.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include <iostream>

// RetinalCamera Definitions
RetinalCamera::RetinalCamera(const AnimatedTransform &cam2world, const float screenWindow[4],
					float sopen, float sclose, Film *film, float pupilradius, float focald, Point nodalP,
					const float *zVals, const Normal *Normals ) 

	: ProjectiveCamera(cam2world, Orthographic(0., 1.), screenWindow, sopen, sclose, 0.0, focald, film)
{
	// compute differential changes in origin for retinal camera rays
	dxCamera = RasterToCamera(Vector(1, 0, 0));
	dyCamera = RasterToCamera(Vector(0, 1, 0));

    numCones    = 131072;
    xResolution = 512;
    yResolution = 256;
    NodalPoint  = nodalP;
    PupilRadius = pupilradius;

    // set the z values and normals of the retina
    if (zVals) 
    {
        zValues = new float[numCones];
        memcpy(zValues, zVals, numCones*sizeof(float));
    }
    if (Normals) 
    {
        rNormals = new Normal[numCones];
        memcpy(rNormals, Normals, numCones*sizeof(Normal));
    }
}

RetinalCamera::~RetinalCamera()
{
    delete[] zValues;
    // delete[] rNormals; NEED TO CHANGE THIS LATER
}


float RetinalCamera::GenerateRay( const CameraSample &sample, Ray * ray) const
{
	// The new Origin =  (sample.imageX, sample.imageY, zVals[sample.imageX + sample.imageY * rowL])
	// the time should be set to 0
    int row = int( sample.imageY );
    int col = int( sample.imageX );
    float z = zValues[ row * xResolution + col ];

    Point Pras(sample.imageX, sample.imageY, z);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);

    Vector Rdir(NodalPoint - Pras);
    Vector Cdir;
    RasterToCamera(Rdir, &Cdir);

    *ray = RayDifferential(Pcamera, Normalize(Cdir), 0., INFINITY);
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
    int row = int( sample.imageY );
    int col = int( sample.imageX );
    float z = zValues[ row * xResolution + col ];

    Point Pras(sample.imageX, sample.imageY, z);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);

    Vector Rdir(NodalPoint - Pras);
    Vector Cdir;
    RasterToCamera(Rdir, &Cdir);

    *ray = RayDifferential(Pcamera, Normalize(Cdir), 0., INFINITY);

    ray->time = sample.time;

    CameraSample sshift = sample;

    Point rxNew(++sshift.imageX, sshift.imageY, zValues[ row * xResolution + col + 1 ] );
    Point ryNew(--sshift.imageX, ++sshift.imageY, zValues[ (row + 1) * xResolution + col ] );
    Vector rxDir = NodalPoint - rxNew;
    Vector ryDir = NodalPoint - ryNew;
    
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
    float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        swap(shutterclose, shutteropen);
    }

    Point nodalpoint  = params.FindOnePoint("nodalpoint", Point(256.0, 128.0, 17.2) );
    float pupilradius = params.FindOneFloat("pupilradius", 0.f);
    float focaldistance = params.FindOneFloat("focaldistance", 1e30f);
    float frame = params.FindOneFloat("frameaspectratio", float(film->xResolution)/float(film->yResolution));

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

    const float *sw = params.FindFloat("screenwindow", &swi);
    const float *zVals = params.FindFloat("zvalues", &zvi);
    const Normal *Normals = params.FindNormal("normals", &rni);

    if (sw && swi == 4)
    {
        memcpy( screen, sw, 4*sizeof(float) );
    }

    
    return new RetinalCamera(cam2world, screen, shutteropen, shutterclose,
        film, pupilradius, focaldistance, nodalpoint, zVals, Normals);
}
