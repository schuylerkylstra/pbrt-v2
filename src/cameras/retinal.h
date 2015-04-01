

// Written by: Schuyler Kylstra
// Comp 992
// Spring 2015
// Feel free to use/alter this code to your whim

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_RETINAL_H
#define PBRT_CAMERAS_RETINAL_H


#include "pbrt.h"
#include "camera.h"
#include "film.h"


class RetinalCamera : public ProjectiveCamera
{
public:
    
    
    RetinalCamera(	const AnimatedTransform &cam2world, const float screenWindow[4],
    				float sopen, float sclose, Film *film, float lensr, float focald, Point nodalP, 
    				const float *zVals, const Normal *Normals );
    ~RetinalCamera();
    float GenerateRay( const CameraSample &sample, Ray * ) const;
    float GenerateRayDifferential( const CameraSample &sample, RayDifferential *ray ) const;
    
    
private:
	// Retinal Camera private data
    float   *zValues; // 1st collumn = z[0], z[256], z[512], etc.
    Normal  *rNormals; 
    int     numCones, xResolution, yResolution;
    Vector  dxCamera, dyCamera;
    Point   NodalPoint;
    float   PupilRadius;
    
};


RetinalCamera* CreateRetinalCamera(const ParamSet &params, const AnimatedTransform &cam2world, Film *film);



#endif