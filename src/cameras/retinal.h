

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
#include <iostream>

struct Lens
{
    Lens(float ri = 1.33, float lr = 5.5, Point o = Point(256.0,128.0, 0.0) )
    {
        RefractiveIndex = ri;
        LensRadius = lr;
        Center = o;
    }
    ~Lens(){}

    float RefractiveIndex;
    float LensRadius;
    Point Center;


    void ApplySnellLaw( Ray * r ) const
    {
        // This will always yeild the farther point but that may not be a good point
        Point Impact;
        if ( Intersection(r, &Impact) )
        {
            // determine the normal
            Vector n(Impact - Center);
            n = Normalize(n);

            // normalize incident direction
            Vector dir = r->d;
            dir = Normalize(dir);

            // calculate transmition direction
            float dot = Dot(dir, n);
            // std::cout << "DOT: " << dot << "\n";

            Vector t = dir * RefractiveIndex + n * ( dot * RefractiveIndex - sqrt(1 - RefractiveIndex * RefractiveIndex * ( 1 - dot * dot) ) );
            // std::cout << "T before: " << t.x << ", " << t.y << ", " << t.z << "\n";
            // std::cout << "Radical: " << 1 - dot * dot * RefractiveIndex * RefractiveIndex << "\n";
            // t = t + n * ( dot * RefractiveIndex - sqrt(1 - RefractiveIndex * RefractiveIndex * ( 1 - dot * dot) ) );
            // std::cout << "T after: " << t.x << ", " << t.y << ", " << t.z << "\n";

            // update the ray
            r->o = Impact;
            r->d = t;
        }
    }

    float quadEQ( float a, float b, float c ) const
    {
        float rad = b*b - 4*a*c;
        if( rad < 0.f)
            return -1;
        rad = sqrt(rad);
        float s1 = (-b + rad)/(2*a);
        float s2 = (-b - rad)/(2*a);
        if( s1 < 0.f && s2 < 0.f)
            return -1;
        if(s1 >= s2)
            return s1;
        return s2;
    }

    bool Intersection(Ray * r, Point * impact) const
    {
        Vector offset(r->o - Center);
        // Dot(r->d,r->d), 2*(Dot(r->d, offset)), Dot(offset,offset) - LensRadius*LensRadius
        float t = quadEQ( Dot(r->d,r->d), 2*(Dot(r->d, offset)), Dot(offset,offset) - LensRadius*LensRadius );
        if ( t > 0.f)
        {
            *impact = (r->o) + (r->d)*t;
            return true;
        }
        return false;
    }


};

class RetinalCamera : public ProjectiveCamera
{
public:
    
    
    RetinalCamera(	const AnimatedTransform &cam2world, const float screenWindow[4],
    				float sopen, float sclose, Film *film, float lensr, float focald, Point nodalP, 
    				const float *zVals, const Normal *Normals , float rIndex);
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
    Lens    Cornea;
    bool    zFlag, nFlag;
    
};


RetinalCamera* CreateRetinalCamera(const ParamSet &params, const AnimatedTransform &cam2world, Film *film);



#endif