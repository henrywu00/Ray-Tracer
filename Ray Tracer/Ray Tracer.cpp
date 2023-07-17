// Ray Tracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _CRT_SECURE_NO_WARNINGS

#define _USE_MATH_DEFINES

#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>

#include "glm/glm.hpp"
#include "glm/ext.hpp"
#include "lodepng.h"
#include "lodepng.cpp"

using namespace std;

float factor = 0.8;

class Camera {
public:
    glm::vec3 lookAt;
    glm::vec3 up;
    glm::vec3 eye;
    int resx, resy, fov;

    Camera(int x, int y) {
        lookAt = glm::vec3(0, 0, 0);
        up = glm::vec3(0, 1, 0);
        eye = glm::vec3(0, 0, -150);
        resx = x;
        resy = y;
        fov = 90;
    }
};

class Ray {
public:
    glm::vec3 origin;
    glm::vec3 dir;

    Ray(glm::vec3 d) {
        origin = glm::vec3(0);
        dir = d;
    }

    Ray(glm::vec3 o, glm::vec3 d) {
        origin = o;
        dir = d;
    }

    void printRay() {
        cout << "origin: " << glm::to_string(origin) << ", dir: " << glm::to_string(dir) << endl;
    }
};

class Light {
public:
    glm::vec3 pos;
    glm::vec3 diff;
    glm::vec3 spec;

    Light(glm::vec3 p, glm::vec3 d, glm::vec3 s) {
        pos = p;
        diff = d;
        spec = s;
    }

    void printLight() {
        cout << glm::to_string(pos) << endl;
        cout << glm::to_string(diff) << endl;
        cout << glm::to_string(spec) << endl << endl;
    }

}; 

class Sphere{
public:
    glm::vec3 pos;
    glm::vec3 diff;
    glm::vec3 spec;
    double shininess;
    double radius;
    bool opaque;
    double index;

    Sphere() {
        glm::vec3 i(0);
        pos = diff = spec = i;
        radius = 1;
        shininess = 0;
        opaque = false;
        index = 0;
    }

    Sphere(glm::vec3 p, double r, glm::vec3 d, glm::vec3 s, double shine, bool op, double idx) {
        pos = p;
        radius = r;
        diff = d;
        spec = s;
        shininess = shine;
        opaque = op;
        index = idx;
    }



    void printSphere(){
        cout << glm::to_string(pos) << endl;
        cout << radius << endl;
        cout << glm::to_string(diff) << endl;
        cout << glm::to_string(spec) << endl;
        cout << shininess << endl;
        cout << opaque << endl;
        cout << index << endl << endl;
    }

    pair<float, float> intersect(Ray ray) {
        double b = 2.0f * (ray.dir.x * (ray.origin.x - pos.x) + ray.dir.y * (ray.origin.y - pos.y) + ray.dir.z * (ray.origin.z - pos.z));
        double c = pow((ray.origin.x - pos.x), 2) + pow((ray.origin.y - pos.y), 2) + pow((ray.origin.z - pos.z), 2) - pow(radius, 2);

        double d = pow(b, 2) - 4.0f * c;

        if (d < 0) return make_pair(-1, -1);

        float t1 = (-b - sqrt(pow(b, 2) - 4 * c)) / 2.0f;
        float t2 = (-b + sqrt(pow(b, 2) - 4 * c)) / 2.0f;

        return make_pair(t1, t2);
    }
    
};

class Quad {
public:
    glm::vec3 a;
    glm::vec3 b;
    glm::vec3 c;
    glm::vec3 diff;
    glm::vec3 spec;
    double shininess;

    Quad() {
        glm::vec3 i(0);
        a = b = c = diff = spec = i;
        shininess = 0;
    }

    Quad(glm::vec3 ap, glm::vec3 bp, glm::vec3 cp, glm::vec3 d, glm::vec3 s, double sh) {
        a = ap;
        b = bp;
        c = cp;
        diff = d;
        spec = s;
        shininess = sh;
    }

    float intersect(Ray ray) {
        glm::vec3 p = b - a;
        glm::vec3 q = c - a;
        glm::vec3 tmp1 = glm::cross(ray.dir, q);
        float dot1 = glm::dot(tmp1, p);
        if (dot1 > -0.0001 && dot1 < 0.0001) return 0;

        float f = 1.0f / dot1;
        glm::vec3 s = ray.origin - a;
        float u = f * glm::dot(s, tmp1);
        if (u < 0.0f || u > 1.0f) return 0;

        glm::vec3 tmp2 = glm::cross(s, p);
        float v = f * glm::dot(ray.dir, tmp2);
        if (v < 0.0f || v > 1.0f) return 0;

        float t = f * glm::dot(q, tmp2);

        return t;
    }

};

class Scene {
public:
    int depth;
    glm::vec3 background;
    int resx;
    int resy;
    vector<Sphere> spheres;
    vector<Quad> quads;
    vector<Light> lights;


    Scene(string path) {
        ifstream infile(path);

        int x, y, z, r;
        double a, b, c, shine;

        infile >> a >> b >> c;
        background = glm::vec3(a, b, c);

        infile >> depth;

        int num;

        bool refrac;

        infile >> num;

        for (int i = 0; i < num; i++) {
            infile >> x >> y >> z;
            glm::vec3 pos = glm::vec3(x, y, z);

            infile >> a >> b >> c;
            glm::vec3 diff = glm::vec3(a, b, c);

            infile >> a >> b >> c;
            glm::vec3 spec = glm::vec3(a, b, c);

            lights.push_back(Light(pos, diff, spec));
        }

        infile >> num;

        for (int i = 0; i < num; i++) {
            infile >> x >> y >> z;
            glm::vec3 pos = glm::vec3(x, y, z);

            infile >> r;
            
            infile >> a >> b >> c;
            glm::vec3 diff = glm::vec3(a, b, c);

            infile >> a >> b >> c;
            glm::vec3 spec = glm::vec3(a, b, c);

            infile >> shine;

            double idx;

            infile >> refrac >> idx;

            spheres.push_back(Sphere(pos, r, diff, spec, shine, refrac, idx));
        }

        infile >> num;

        for (int i = 0; i < num; i++) {
            infile >> x >> y >> z;
            glm::vec3 ap = glm::vec3(x, y, z);

            infile >> x >> y >> z;
            glm::vec3 bp = glm::vec3(x, y, z);

            infile >> x >> y >> z;
            glm::vec3 cp = glm::vec3(x, y, z);

            infile >> a >> b >> c;
            glm::vec3 diff = glm::vec3(a, b, c);

            infile >> a >> b >> c;
            glm::vec3 spec = glm::vec3(a, b, c);

            infile >> shine;

            quads.push_back(Quad(ap, bp, cp, diff, spec, shine));
        }

        infile >> resx >> resy;

    }

    void printScene() {
        cout << glm::to_string(background) << endl;
        cout << depth << endl;
        cout << resx << " " << resy << endl;
        for (Light light : lights) {
            light.printLight();
        }
        for (Sphere sphere : spheres) {
            sphere.printSphere();
        }
    }
};


Ray calculateRay(int i, int j, Camera camera) {
    glm::vec3 l = (camera.lookAt - camera.eye) / glm::length(camera.lookAt - camera.eye);
    glm::vec3 v = glm::cross(l, camera.up) / glm::length(glm::cross(l, camera.up));
    glm::vec3 u = glm::cross(v, l);

    float a = camera.resx / camera.resy;
    float d = 1 / tan((camera.fov / 2) * (M_PI / 180.0f));

    glm::vec3 ll = camera.eye + l * d - v * a + u;

    glm::vec3 p = ll + (2 * a * v * (float)i) / (float)camera.resx - 2.0f * u * (float) j / (float) camera.resy;
    glm::vec3 dd = (p - camera.eye) / glm::length(p - camera.eye);
    return Ray(camera.eye, dd);
}

pair<glm::vec3, Sphere> firstSphereIntersection(Ray ray, Scene scene) {
    glm::vec3 p(0);
    if (!scene.spheres.empty()) {
        Sphere s = scene.spheres[0];
        double min = INT_MAX;
        for (Sphere sphere : scene.spheres) {
            pair<float, float> curr = sphere.intersect(ray);

            if (curr.first > 0) {

                glm::vec3 point = ray.origin + ray.dir * curr.first;

                if (glm::length(point) < min) {
                    min = glm::length(point);
                    p = point;
                    s = sphere;
                }
            }
        }
        return make_pair(p, s);
    }
    else return make_pair(p, Sphere());
    
}

pair<glm::vec3, Quad> firstQuadIntersection(Ray ray, Scene scene) {
    glm::vec3 p(0);
    if (!scene.quads.empty()) {
        Quad q = scene.quads[0];
        double min = INT_MAX;
        for (Quad quad : scene.quads) {
            float curr = quad.intersect(ray);

            if (curr > 0) {
                glm::vec3 point = ray.origin + ray.dir * curr;

                if (glm::length(point) < min) {
                    min = glm::length(point);
                    p = point;
                    q = quad;
                }
            }
        }
        return make_pair(p, q);
    }
    else return make_pair(p, Quad());
    
    
}

vector<Light> shadowRaysS(pair<glm::vec3, Sphere> point, Scene scene) {
    //cout << "first intersection: " << glm::to_string(point.first) << endl;
    glm::vec3 norm = glm::normalize(point.first - point.second.pos);

    pair<float, float> curr = make_pair(-1, -1);;
    vector<Light> visible;
    for (Light light : scene.lights) {
        glm::vec3 lightDir = glm::normalize(light.pos - (point.first + norm * 0.1f));
        
        Ray ray(point.first + norm * 0.1f, lightDir);
        double dist = glm::length(light.pos - (point.first + norm * 0.1f));
        for (Sphere sphere : scene.spheres) {

            curr = sphere.intersect(ray);
            //cout << "t1: " << curr.first << ", t2: " << curr.second << endl;
            if (curr.first > 0 && curr.second > 0) {
                //cout << "intersected!" << endl;
                if (curr.first < dist || curr.second < dist) {
                    break;                  
                }
            }
        }

        float currQ = 0;

        for (Quad quad : scene.quads) {

            currQ = quad.intersect(ray);
            //cout << "t1: " << curr.first << ", t2: " << curr.second << endl;
            if (currQ > 0) {
                //cout << "intersected!" << endl;
                if (currQ < dist) break;
            }
        }

        if ((curr.first <= 0 && curr.second <= 0 && currQ <= 0) || (curr.first > dist && curr.second > dist || currQ > dist)) {
            visible.push_back(light);
        }
    }
    return visible;
}

vector<Light> shadowRaysQ(pair<glm::vec3, Quad> point, Scene scene) {

    glm::vec3 norm = glm::normalize(glm::cross(point.second.b - point.second.a, point.second.c - point.second.a));

    pair<float, float> curr = make_pair(-1, -1);
    vector<Light> visible;
    for (Light light : scene.lights) {
        glm::vec3 lightDir = glm::normalize(light.pos - (point.first + norm * 0.001f));

        Ray ray(point.first + norm * 0.001f, lightDir);
        double dist = glm::length(light.pos - (point.first + norm * 0.001f));
        for (Sphere sphere : scene.spheres) {

            curr = sphere.intersect(ray);

            if (curr.first > 0 && curr.second > 0) {

                if (curr.first < dist || curr.second < dist) {
                    break;
                }
            }
        }

        float currQ = 0;

        for (Quad quad : scene.quads) {
            currQ = quad.intersect(ray);

            if (currQ > 0) {

                if (currQ < dist) break;
                
            }
        }

        if ((curr.first <= 0 && curr.second <= 0 && currQ <= 0) || (curr.first > dist && curr.second > dist || currQ > dist)) {
            visible.push_back(light);
        }
    }
    return visible;
}

glm::vec3 phongS(pair<glm::vec3, Sphere> intersection, Light light, Ray ray) {
    glm::vec3 l = glm::normalize(light.pos - intersection.first);
    glm::vec3 n = glm::normalize(intersection.first - intersection.second.pos);
    glm::vec3 dmult = glm::vec3(intersection.second.diff.x * light.diff.x, intersection.second.diff.y * light.diff.y, intersection.second.diff.z * light.diff.z);

    glm::vec3 smult = glm::vec3(intersection.second.spec.x * light.spec.x, intersection.second.spec.y * light.spec.y, intersection.second.spec.z * light.spec.z);

    glm::vec3 v = glm::normalize(intersection.first - ray.origin);
    glm::vec3 r = glm::reflect(l, n);
    float sDot = max(glm::dot(l, n), 0.0f);

    glm::vec3 spec = glm::vec3(0);

    if (sDot > 0) {
        spec = smult * (float)pow(max(glm::dot(r, v), 0.0f), intersection.second.shininess);
    }

    glm::vec3 diffuse = dmult * sDot;
    return diffuse + spec;
}

glm::vec3 phongQ(pair<glm::vec3, Quad> intersection, Light light, Ray ray) {
    glm::vec3 l = glm::normalize(light.pos - intersection.first);
    glm::vec3 n = glm::normalize(glm::cross(intersection.second.b - intersection.second.a, intersection.second.c - intersection.second.a));
    glm::vec3 dmult = glm::vec3(intersection.second.diff.x * light.diff.x, intersection.second.diff.y * light.diff.y, intersection.second.diff.z * light.diff.z);

    glm::vec3 smult = glm::vec3(intersection.second.spec.x * light.spec.x, intersection.second.spec.y * light.spec.y, intersection.second.spec.z * light.spec.z);

    glm::vec3 v = glm::normalize(intersection.first - ray.origin);
    glm::vec3 r = glm::reflect(l, n);
    float sDot = max(glm::dot(l, n), 0.0f);

    glm::vec3 spec = glm::vec3(0);

    if (sDot > 0) {
        spec = smult * (float)pow(max(glm::dot(r, v), 0.0f), intersection.second.shininess);
    }

    glm::vec3 diffuse = dmult * sDot;
    return diffuse + spec;
}

Ray reflectedS(Ray ray, glm::vec3 point, Sphere sphere) {
    glm::vec3 n = glm::normalize(point - sphere.pos);
    glm::vec3 r = glm::reflect(ray.dir, n);
    return Ray(point, r);
}

Ray reflectedQ(Ray ray, glm::vec3 point, Quad quad) {
    glm::vec3 n = glm::normalize(glm::cross(quad.b - quad.a, quad.c - quad.a));
    glm::vec3 r = glm::reflect(ray.dir, n);
    return Ray(point + n * 0.1f, r);
}

Ray refracted(Ray ray, glm::vec3 point, Sphere sphere, float inside) {
    glm::vec3 n = glm::normalize(point - sphere.pos) * inside;
    double n1 = 1.0;
    double n2 = sphere.index;
    glm::vec3 v = (float)(n1 / n2) * ((float)(sqrt(pow(glm::dot(n, ray.dir), 2) + pow(n2 / n1, 2) - 1) - glm::dot(n, ray.dir)) * n + ray.dir);
    return Ray(point, glm::normalize(v));
}

glm::vec3 traceRay(Ray ray, int depth, Scene scene, float inside) {
    glm::vec3 color = scene.background;

    pair<glm::vec3, Sphere> interS = firstSphereIntersection(ray, scene);

    pair<glm::vec3, Quad> interQ = firstQuadIntersection(ray, scene);

    float distS = glm::length(interS.first - ray.origin);
    float distQ = glm::length(interQ.first - ray.origin);

    if (interS.first != glm::vec3(0)) {
        if ((interQ.first != glm::vec3(0) && (distS < distQ)) || (interQ.first == glm::vec3(0))) {
            color += interS.second.diff * 0.05f;

            vector<Light> contributingLights = shadowRaysS(interS, scene);
            if (depth <= 0) return color;

            for (Light light : contributingLights) {
                color += phongS(interS, light, ray) * factor;
            }

            color += traceRay(reflectedS(ray, interS.first, interS.second), depth - 1, scene, inside) * (1 - factor);
        }
    }
    else if (interQ.first != glm::vec3(0)) {
        color += interQ.second.diff * 0.05f;
        vector<Light> contributingLights = shadowRaysQ(interQ, scene);
        if (depth <= 0) return color;

        for (Light light : contributingLights) {
            color += phongQ(interQ, light, ray) * factor;
        }

        //color += traceRay(reflectedQ(ray, interQ.first, interQ.second), depth - 1, scene, inside) * (1 - factor);
    }

    return color;

    
}


void encodeImage(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
    unsigned error = lodepng::encode(filename, image, width, height);

    if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

int main()
{
    Scene scene = Scene("./cornell.txt");
    scene.printScene();


    string name = "test-";
    string end = ".png";

    unsigned width = scene.resx, height = scene.resy;
    std::vector<unsigned char> image;

    Camera camera(width, height);

    image.resize(width * height * 4);


    int x0 = 0;
    int y0 = 0;
    int z0 = 100;
    int r = 80;

    glm::vec3 l = glm::vec3(-1, 1, 1) / glm::length(glm::vec3(-1, 1, 1));
    glm::vec3 up = glm::vec3(0, 1, 0);
    glm::vec3 v = glm::cross(l, up) / glm::length(glm::cross(l,up));
    glm::vec3 u = glm::cross(v, l) / glm::length(glm::cross(v, l));

    string filename = "";

    //This is for rendering one image in the scene
    #pragma omp parallel for
    for (unsigned x = 0; x < width; x++) {
        for (unsigned y = 0; y < height; y++) {

            Ray ray = calculateRay(x, y, camera);

            glm::vec3 color = traceRay(ray, scene.depth, scene, 1);

            image[4 * width * y + 4 * x + 0] = min(color.x * 255, 255.0f);
            image[4 * width * y + 4 * x + 1] = min(color.y * 255, 255.0f);
            image[4 * width * y + 4 * x + 2] = min(color.z * 255, 255.0f);

            image[4 * width * y + 4 * x + 3] = 255;
        }
    }
    filename = name + end;
    encodeImage(filename.c_str(), image, width, height);


    //This commented block is used to create 360 images to make an animation. In this case I revolved a light source around the scene.

    /*for (int t = 0; t < 360; t++) {
        int x = x0 + r * cos(t * M_PI / 180) * v.x + r * sin(t * M_PI / 180) * u.x;
        int y = y0 + r * cos(t * M_PI / 180) * v.y + r * sin(t * M_PI / 180) * u.y;
        int z = z0 + r * cos(t * M_PI / 180) * v.z + r * sin(t * M_PI / 180) * u.z;
        scene.lights[0].pos = glm::vec3(x, y, z);

        #pragma omp parallel for
        for (unsigned x = 0; x < width; x++) {
            for (unsigned y = 0; y < height; y++) {

                Ray ray = calculateRay(x, y, camera);

                glm::vec3 color = traceRay(ray, scene.depth, scene, 1);

                image[4 * width * y + 4 * x + 0] = min(color.x * 255, 255.0f);
                image[4 * width * y + 4 * x + 1] = min(color.y * 255, 255.0f);
                image[4 * width * y + 4 * x + 2] = min(color.z * 255, 255.0f);

                image[4 * width * y + 4 * x + 3] = 255;
            }
        }
        filename = name + to_string(t) + end;
        encodeImage(filename.c_str(), image, width, height);
        cout << "rendered frame " << t << endl;
    }*/

    

}