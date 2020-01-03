#pragma once

#include "ShaderFlat.h"
#include "Scene.h"

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

        /**
         * For Reflection
        */
        Vec3f resReflect(0, 0, 0);
        bool onReflection = false;

        // now we use the reflection vector to create a new ray at the point of intersection:
        if(!isOpaque) {
        // if(!isOpaque && 0==1) {
            onReflection = true;

            Ray rayReflect;
            rayReflect.org = ray.org + ray.t * ray.dir;
            rayReflect.dir = reflect;
            rayReflect.t = std::numeric_limits<float>::infinity();
            rayReflect.reflectDepth = ray.reflectDepth + 1;

            if(ray.reflectDepth <= MAX_REFLECT_DEPTH) {
                resReflect = m_scene.RayTrace(rayReflect);
            }else {
                onReflection = false;
            }

            if(rayReflect.t == std::numeric_limits<float>::infinity())
                onReflection = false;
        }

        if(onReflection && ray.reflectDepth != 0) {
            return resReflect;
        }



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
                            if(shadow.hit->getShader().get()->getIsOpaque()) 
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

        if(onReflection && ray.reflectDepth == 0) {
            res = resReflect + 0.2 * res;
        }

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
    const int MAX_REFLECT_DEPTH = 2;
    const int MAX_REFRACT_DEPTH = 1;
};
