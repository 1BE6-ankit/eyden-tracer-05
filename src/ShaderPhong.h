#pragma once

#include "ShaderFlat.h"
#include "Scene.h"

// #define REFLECT_OFF
// #define REFRACT_OFF

const int nAreaSamples = 1;

class CShaderPhong : public CShaderFlat
{
public:
	/**
	 * @brief Constructor
	 * @param scene Reference to the scene
	 * @param color The color of the object
	 * @param ka The ambient coefficient
	 * @param kd The diffuse reflection coefficients
	 * @param ks The specular refelection coefficients
	 * @param ke The shininess exponent
	 */
	CShaderPhong(CScene& scene, Vec3f color, float ka, float kd, float ks, float ke, bool isOpaque = true)
		: CShaderFlat(color, isOpaque)
		, m_scene(scene)
		, m_ka(ka)
		, m_kd(kd)
		, m_ks(ks)
		, m_ke(ke)
	{}
	virtual ~CShaderPhong(void) = default;

	virtual Vec3f Shade(const Ray& ray) const override
	{
		// get shading normal
		Vec3f normal = ray.hit->getNormal(ray);

		// turn normal to front
		if (normal.dot(ray.dir) > 0)
			normal = -normal;

		// calculate reflection vector
		Vec3f reflect = normalize(ray.dir - 2 * normal.dot(ray.dir) * normal);

		// ambient term
		Vec3f ambientIntensity(1,1,1);

		Vec3f color = CShaderFlat::Shade();
		Vec3f ambientColor = m_ka * color;
		Vec3f res = ambientColor.mul(ambientIntensity);

		// shadow ray (up to now only for the light direction)
		Ray shadow;
		shadow.org = ray.org + ray.t * ray.dir;

		// iterate over all light sources
		for (auto pLight : m_scene.m_vpLights)
			for(int s = 0; s < nAreaSamples; s++) {	// TODO: make the sampling to depend on the light source type
				// get direction to light, and intensity
				std::optional<Vec3f> lightIntensity = pLight->Illuminate(shadow);
				if (lightIntensity) {
					// diffuse term
					float cosLightNormal = shadow.dir.dot(normal);
					if (cosLightNormal > 0) {
						if (m_scene.Occluded(shadow)) {
                            continue;
                        }

						Vec3f diffuseColor = m_kd * color;
						res += (diffuseColor * cosLightNormal).mul(lightIntensity.value());
					}

					// specular term
					float cosLightReflect = shadow.dir.dot(reflect);
					if (cosLightReflect > 0) {
						Vec3f specularColor = m_ks * RGB(1, 1, 1); // white highlight;
						res += (specularColor * powf(cosLightReflect, m_ke)).mul(lightIntensity.value());
					}
				}
			}

		if (nAreaSamples > 1)
			res /= nAreaSamples;

        /**
         * For Refraction
        */
        float nAir = 1, nGlass = 1.5;
        Vec3f resRefract(0, 0, 0);
        bool onRefraction = false;
        bool onIntReflection = false; // total internal reflection flag

#ifndef REFRACT_OFF
        if(!isOpaque && ray.refractDepth <= MAX_REFRACT_DEPTH) {
            onRefraction = true;
            Vec3f refrNormal = normal;
            float nDotI = refrNormal.dot(ray.dir); 

            // since normal was negated in the beginning, nDoI is always < 0
            nDotI = -nDotI;

            if(ray.refractDepth % 2 == 1) 
                std::swap(nAir, nGlass);

            float n = nAir / nGlass;

            Ray refractRay;
            refractRay.org = ray.org + ray.t * ray.dir;
            refractRay.t = std::numeric_limits<float>::infinity();

            // part of formula to calculate refracted ray direction
            float sqrtVal = sqrt(1 - pow(n, 2) * (1 - pow(nDotI, 2))); 

            if(isnan(sqrtVal) ) {
                // total internal reflection
                onIntReflection = true;
                refractRay.dir = reflect;
                refractRay.refractDepth = ray.refractDepth + 2;
                refractRay.reflectDepth = ray.reflectDepth + 1;
            } else {
                // use the formula to calculate direction vector of the refraction ray
                refractRay.dir = normalize(n * (ray.dir + refrNormal * nDotI) - refrNormal * sqrtVal);
                refractRay.refractDepth = ray.refractDepth + 1;
                refractRay.reflectDepth = ray.reflectDepth;
            }

            resRefract = m_scene.RayTrace(refractRay);
        }

        #ifdef REFLECT_OFF
            // For Problem 5.2
            if(onRefraction && ray.refractDepth != 0) {
                return res * (1 - c_transmit) + resRefract * c_transmit;
            }
        #endif
#endif

        /**
         * For Reflection
        */
        Vec3f resReflect(0, 0, 0);
        bool onReflection = false;
#ifndef REFLECT_OFF
        if(!isOpaque && ray.reflectDepth <= MAX_REFLECT_DEPTH && !onIntReflection) {
            onReflection = true;

            Ray rayReflect;
            rayReflect.org = ray.org + ray.t * ray.dir;
            rayReflect.dir = reflect;
            rayReflect.t = std::numeric_limits<float>::infinity();
            rayReflect.reflectDepth = ray.reflectDepth + 1;
            rayReflect.refractDepth = ray.refractDepth;

            resReflect =  m_scene.RayTrace(rayReflect);
        }

        #ifdef REFRACT_OFF
            // For Problem 5.1
            if(onReflection && ray.reflectDepth != 0) {
                return resReflect;
            }
        #endif
#endif

        Vec3f resTemp(0, 0, 0);

#if !defined(REFLECT_OFF) && !defined(REFRACT_OFF)
        // Reflection and Refraction at depth > 0. 
        if(onReflection || onRefraction) {
            if(ray.reflectDepth > 0 || ray.refractDepth > 0)
                return (1 - c_transmit) * res +  c_reflect * resReflect + c_transmit * resRefract;
        }
#endif

#ifndef REFLECT_OFF
        // handle reflection shade at recursion base
        if(onReflection && ray.reflectDepth == 0) {
            if(resReflect.val[0] == 0 && resReflect.val[1] == 0 && resReflect.val[2] == 0) {
                resTemp += res * (1 - c_reflect);
            } else {
                #if defined(REFRACT_OFF)
                    resTemp += resReflect * (1 - c_reflect);
                #else
                    resTemp += resReflect * c_reflect;
                #endif
            }
        }
#endif

#ifndef REFRACT_OFF
        // handle refraction shade at recursion base
        if(onRefraction && ray.refractDepth == 0) {
            if(resRefract.val[0] == 0 && resRefract.val[1] == 0 && resRefract.val[2] == 0) {
                resTemp += res * c_transmit;
            } else {
                resTemp += resRefract * c_transmit + res * (1 - c_transmit);
            }
        } 
#endif

        if(onReflection || onRefraction)
            res = resTemp;

		for (int i = 0; i < 3; i++)
			if (res.val[i] > 1) res.val[i] = 1;
		
		return res;
	}

	
private:
	CScene& m_scene;
	float 	m_ka;    ///< ambient coefficient
	float 	m_kd;    ///< diffuse reflection coefficients
	float 	m_ks;    ///< specular refelection coefficients
	float 	m_ke;    ///< shininess exponent

    const float c_reflect = 0.2; // reflection coefficient
    const float c_transmit = 0.8; // transmission coefficient
    const int MAX_REFLECT_DEPTH = 3;
    const int MAX_REFRACT_DEPTH = 5;
};
