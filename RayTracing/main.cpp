#include "rtweekend.h"
#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "vec3.h"
#include "material.h"
#include <iostream>
#include <fstream>

//Path: G:\Documents\KTH\DD1395\RayTracing\Release

double lengthS(point3 p);
double length_squaredS(point3 p);

double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius * radius;
    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return (-half_b - sqrt(discriminant)) / a;
    }
}


color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth - 1);
        return color(0, 0, 0);
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1, 1, 1) + t * color(0.1, 0.5, 0.5);
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double lengthS(point3 p) {
    return sqrt(length_squaredS(p));
}

double length_squaredS(point3 p){
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
}

hittable_list generateRandom() {
    hittable_list world;
    
    float length = 4;

    for (float x = -length; x <= length; x+=0.5) {
        for (float z = -length; z <= length; z+=0.5) {

            float y = 2 * x * x + 2 * z * z - 4;
            auto choice = (int)ceil(fRand(0, 3));
            float size = 0.2;//fmax(t / 12,0.5);
            std::cerr << x << " " << y << " " << z << " " << std::endl;
            if (choice == 1) {
                //regular

                auto material = make_shared<lambertian>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)));
                world.add(make_shared<sphere>(point3(x, y, z), size, material));
            }
            else if (choice == 2)
            {

                auto material = make_shared<metal>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)), fRand(0, 1));
                world.add(make_shared<sphere>(point3(x, y, z), size, material));
            }
            else if (choice == 3)
            {
                auto material = make_shared<dielectric>(fRand(-0.5, 0.5));
                world.add(make_shared<sphere>(point3(x, y, z), size, material));
            }
        }

    }



    return world;
}


hittable_list generateSpiral() {
    hittable_list world;
    int n = 100;
    //int t = 0;

    for (float t = 0; t <= 50; t = t + fRand(0.1,0.2)) {
        std::cerr << "\rHey! " << t << ' ' << std::endl;
        auto x = t * cos(3 * t);
        auto z = t * sin(3 * t);
        auto y = t;


        auto choice = (int)ceil(fRand(0, 3));
        float size = .6;//fmax(t / 12,0.5);
        std::cerr << size << std::endl;
        if (choice == 1) {
            //regular

            auto material = make_shared<lambertian>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)));
            world.add(make_shared<sphere>(point3(x, y, z), size, material));
        }
        else if (choice == 2)
        {

            auto material = make_shared<metal>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)), fRand(0, 1));
            world.add(make_shared<sphere>(point3(x, y, z), size, material));
        }
        else if (choice == 3)
        {
            auto material = make_shared<dielectric>(fRand(-0.5, 0.5));
            world.add(make_shared<sphere>(point3(x, y, z), size, material));
        }

    }
        


    return world;
}

hittable_list generateGrid() {
    hittable_list world;
    
    //Want to make a grid of spheres with different materials
    double fraction = 10 / (7*6);
    double counter = -1.0;
    int diec = 0;
    for (int i = 9; i >=-9; i-=3)
    {
        for (int j = -6; j <= 9; j += 3)
        {
            int choice = (int)counter % 3 +1;//(int)ceil(fRand(0, 2));
            
            auto size = 1;
            if (choice == 1) {
                //regular

                auto material = make_shared<lambertian>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)));
                world.add(make_shared<sphere>(point3(i,j,0), size, material));
            }
            else if (choice == 2)
            {

                auto material = make_shared<metal>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)), fRand(0, 1));
                world.add(make_shared<sphere>(point3(i, j, 0), size, material));
            }
            else/* if (choice == 3)*/
            {
                auto material = make_shared<dielectric>((double) ((int)diec%6) +0.1);
                std::cerr << "\rHey! " << (double)((int)diec % 6) << std::endl;
                diec++;
                
                world.add(make_shared<sphere>(point3(i, j, 0), size, material));
            }
            
            counter++;
        }
    }

    return world;
}

int main() {

    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 360;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    const int max_depth = 40;

    // Camera
    hittable_list world;

    auto material_ground = make_shared<lambertian>(color(1,1,1));
    auto material_center = make_shared<lambertian>(color(1,0,1));
    auto material_left = make_shared<dielectric>(1.5);
    auto material_right = make_shared<metal>(color(.5,.5,.5), 0.0);
    auto material_1 = make_shared<metal>(color(1,1,1), 0);

    world.add(make_shared<sphere>(point3(0.0, -1001, -1.0), 1000.0, material_ground));
    world.add(make_shared<sphere>(point3(0.0, 1007, -1.0), 1000.0, material_ground));
    //world.add(make_shared<sphere>(point3(-125.0, 0, -125.0), 100.0, material_right));
    //world.add(make_shared<sphere>(point3(0.0, 1200.5, -1.0), 1000.0, material_ground));
    //world.add(make_shared<sphere>(point3(0.0, 5, 0), 3, material_right));
    //world.add(make_shared<sphere>(point3(0.0, 0.0, -1.0), 0.5, material_center));
    //world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
    //world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), -0.45, material_left));
    //world.add(make_shared<sphere>(point3(1.0, 0.0, -1.0), 0.5, material_right));

    //world.add(make_shared<sphere>(point3(0, 1, 0), 2, make_shared<dielectric>(2)));
    //world.add(make_shared<sphere>(point3(2, 2, 3), 3, make_shared<dielectric>(1.01)));
    world.add(make_shared<sphere>(point3(7, 4, 0), 5, make_shared<lambertian>(color(1, 0, 1))));
    /*5, 4, -5), 5*/

    //Generate world;
    //world = generateGrid();
    //world = generateRandom();

    /*
    point3 focusPoint(0, 0, 0);
    int spheres = 200;
    for (size_t i = 0; i < spheres; i++)
    {
        std::cerr << "\Spheres remaining: " << i << ' ' << std::flush;
        //point3 pos(fRand(30, 30), fRand(1, 1), fRand(30, 30));
        //add to focus
        //focusPoint += pos;
        //color c(fRand(0, 1), fRand(0, 1), fRand(0, 1));
        auto choice = ceil(fRand(0, 3));
        auto size = fRand(0.1, 2);
        if (choice <= 1) {
            //regular

            auto material= make_shared<lambertian>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)));
            world.add(make_shared<sphere>(point3(fRand(-20,20), fRand(size/2, 5+size/2), fRand(-20, 20)), size, material));
        }
        else if(choice >= 1 && choice
             < 2)
        {
            std::cerr << "\rHey! " << ' ' << std::endl;
            auto material = make_shared<dielectric>(fRand(-5, 1));
            world.add(make_shared<sphere>(point3(fRand(-20, 20), fRand(size / 2, 5 + size / 2), fRand(-20, 20)), size, material));
            
        }
        else        {
            auto material = make_shared<metal>(color(fRand(0, 1), fRand(0, 1), fRand(0, 1)), fRand(0, 1));
            world.add(make_shared<sphere>(point3(fRand(-20, 20), fRand(size / 2, 5 + size / 2), fRand(-20, 20)), size, material));
        }

    }
    focusPoint /= spheres;*/

    // Camera

    point3 lookfrom(-10,4,0);
    point3 lookat(vec3(2.5,4, 0));
 
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 20;
    auto aperture = 0.2;
    
    double distance_Test = lengthS(lookat -lookfrom);

    camera cam(lookfrom, lookat, vup, 50, aspect_ratio, aperture, distance_Test);

    // Render

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    //File
   

    for (int j = image_height - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);

        }
    }

    std::cerr << "\nDone.\n";
}

