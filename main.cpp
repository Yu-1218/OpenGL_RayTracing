#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray &ray,
               const std::vector<LightSource> &light_sources,
               const Hittable &scene,
               int trace_depth);

Color Shade(const std::vector<LightSource> &light_sources,
            const Hittable &hittable_collection,
            const HitRecord &hit_record,
            int trace_depth)
{
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);

    // Ambient term
    color = hit_record.material.k_a * hit_record.material.ambient;

    // Diffuse and Specular terms
    for (int i = 0; i < light_sources.size(); i++)
    {
        Ray raw_shadow_ray = Ray(hit_record.position, glm::normalize(light_sources[i].position - hit_record.position));
        float EER = 1e-3;
        Ray shadow_ray = Ray(hit_record.position + hit_record.normal * EER, raw_shadow_ray.d);

        if (glm::dot(hit_record.normal, shadow_ray.d) > 0)
        {
            HitRecord tmp;

            Vec lightdir = light_sources[i].position - hit_record.position;
            float lightDistance2 = glm::dot(lightdir, lightdir);
            // is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
            auto shadow_res = hittable_collection.Hit(shadow_ray, &tmp);
            bool inShadow = shadow_res && (tmp.distance * tmp.distance < lightDistance2);

            if (!inShadow)
            {
                Vec reflected_shadow_ray = -2 * glm::dot(shadow_ray.d, hit_record.normal) * hit_record.normal + shadow_ray.d;
                Color diffuse = hit_record.material.k_d * hit_record.material.diffuse * glm::dot(hit_record.normal, shadow_ray.d);
                Color specular = hit_record.material.k_s * hit_record.material.specular * pow(glm::dot(reflected_shadow_ray, -1.f * hit_record.in_direction), hit_record.material.sh);
                for (int i = 0; i < 3; i++)
                {
                    if (diffuse[i] < 0)
                    {
                        diffuse[i] = 0;
                    }
                    if (specular[i] < 0)
                    {
                        specular[i] = 0;
                    }
                }
                color += light_sources[i].intensity * (diffuse + specular);
            }
        }
    }

    // Refelcted Ray
    if (trace_depth < kMaxTraceDepth)
    {
        if (hit_record.material.k_s > 0)
        {
            Ray reflected_ray = Ray(hit_record.position, hit_record.reflection);
            Color r_color = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth + 1);
            color += hit_record.material.k_s * r_color;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        if (color[i] > 1.0)
        {
            color[i] = 1.0;
        }
    }

    return color;
}

Color TraceRay(const Ray &ray,
               const std::vector<LightSource> &light_sources,
               const Hittable &hittable_collection,
               int trace_depth)
{
    // TODO: Add your code here.
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);

    if (hittable_collection.Hit(ray, &record))
    {
        color = Shade(light_sources, hittable_collection, record, trace_depth);
        return color;
    }
    else
    {
        return color;
    }
}

int main()
{
    // TODO: Set your workdir (absolute path) here.
    const std::string work_dir("C:\\Users\\41449\\Desktop\\COMP3271\\A3\\Graphics_PA3_Release\\");

    // Construct scene
    Scene scene(work_dir, "scene\\spheres.toml");
    const Camera &camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f)
            {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f)
                {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++)
            {
                image[idx + i] = (uint8_t)(glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f)
            {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "output.png");
}
