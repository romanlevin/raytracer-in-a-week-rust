extern crate rand;

use std::env;

mod vec3;

use crate::rand::Rng;
use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use vec3::Vec3;

#[derive(Copy, Clone, Debug)]
struct Ray {
    origin: Vec3,
    direction: Vec3,
}

impl Ray {
    fn new(origin: Vec3, direction: Vec3) -> Ray {
        Ray { origin, direction }
    }
    fn point_at_parameter(&self, t: f64) -> Vec3 {
        self.origin + t * self.direction
    }
}

#[derive(Clone, Copy, Debug)]
struct Sphere {
    center: Vec3,
    radius: f64,
}

impl Sphere {
    fn new(center: Vec3, radius: f64) -> Sphere {
        Sphere { center, radius }
    }
}

#[derive(Copy, Clone, Debug)]
struct HitRecord {
    t: f64,
    p: Vec3,
    normal: Vec3,
}

trait Hit {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

impl Hit for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(&ray.direction);
        let b = oc.dot(&ray.direction);
        let c = oc.dot(&oc) - self.radius * self.radius;
        let discriminant = b * b - a * c;
        if discriminant < 0.0 {
            return None;
        }
        let temp = -(b + (b * b - a * c).sqrt()) / a;
        if temp < t_max && temp > t_min {
            let point = ray.point_at_parameter(temp);
            return Some(HitRecord {
                t: temp,
                p: point,
                normal: (point - self.center) / self.radius,
            });
        }
        let temp = (-b + (b * b - a * c).sqrt()) / a;
        if temp < t_max && temp > t_min {
            let point = ray.point_at_parameter(temp);
            return Some(HitRecord {
                t: temp,
                p: point,
                normal: (point - self.center) / self.radius,
            });
        }
        None
    }
}

enum Hittable {
    Sphere(Sphere),
}

impl Hittable {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        match self {
            Hittable::Sphere(sphere) => sphere.hit(ray, t_min, t_max),
        }
    }
}

struct HittableList {
    list: Vec<Hittable>,
}

impl<'a> HittableList {
    fn hit_all(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest_so_far = t_max;
        let mut possible_hit_record: Option<HitRecord> = None;

        self.list.iter().for_each(|hittable| {
            if let Some(record) = hittable.hit(ray, t_min, closest_so_far) {
                closest_so_far = record.t;
                possible_hit_record = Some(record)
            }
        });
        possible_hit_record
    }
}

#[derive(Clone, Copy, Debug)]
struct Camera {
    lower_left_corner: Vec3,
    horizontal: Vec3,
    vertical: Vec3,
    origin: Vec3,
}

impl Camera {
    fn new() -> Camera {
        Camera {
            lower_left_corner: Vec3::new(-2.0, -1.0, -1.0),
            horizontal: Vec3::new(4.0, 0.0, 0.0),
            vertical: Vec3::new(0.0, 2.0, 0.0),
            origin: Vec3::new(0.0, 0.0, 0.0),
        }
    }

    fn get_ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(
            self.origin,
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin,
        )
    }
}

fn random_in_unit_sphere() -> Vec3 {
    let mut rng = rand::thread_rng();
    let mut p: Vec3;
    loop {
        p =
            2.0 * Vec3 {
                x: rng.gen(),
                y: rng.gen(),
                z: rng.gen(),
            } - Vec3::new(1.0, 1.0, 1.0);
        if p.squared_length() >= 1.0 {
            return p;
        }
    }
}

fn color(ray: &Ray, world: &HittableList) -> Vec3 {
    if let Some(hit_record) = world.hit_all(&ray, 0.001, std::f64::MAX) {
        let target = hit_record.p + hit_record.normal + random_in_unit_sphere();
        0.5 * color(
            &Ray {
                origin: hit_record.p,
                direction: target - hit_record.p,
            },
            world,
        )
    } else {
        let unit_direction = ray.direction.unit();
        let t = (unit_direction.y + 1.0) * 0.5;
        (1.0 - t) * Vec3::new(1.0, 1.0, 1.0) + t * Vec3::new(0.5, 0.7, 1.0)
    }
}

#[derive(Debug, Clone, Copy)]
struct ImageParams {
    width: u32,
    height: u32,
    samples: u32,
}

fn parse_args() -> ImageParams {
    let args: Vec<String> = env::args().collect();

    match args.len() {
        1 => panic!("Need to pass some arguments!"),
        3 => {
            let parsed_args: Vec<u32> = args[1..]
                .iter()
                .map(|arg| arg.parse::<u32>().unwrap())
                .collect();
            ImageParams {
                height: parsed_args[0],
                width: parsed_args[0] * 2,
                samples: parsed_args[1],
            }
        }
        _ => panic!("No idea!"),
    }
}

fn main() {
    let image_params = parse_args();

    println!("P3");
    println!("{} {}", image_params.width, image_params.height);
    println!("255");

    let camera = Camera::new();

    let world = &HittableList {
        list: vec![
            Hittable::Sphere(Sphere::new(Vec3::new(0.0, 0.0, -1.0), 0.5)),
            Hittable::Sphere(Sphere::new(Vec3::new(0.0, -100.5, -1.0), 100.0)),
        ],
    };

    let colors: Vec<Vec3> = (0..image_params.height)
        .into_par_iter()
        .rev()
        .flat_map(|j| {
            (0..image_params.width).into_par_iter().map(move |i| {
                let mut rng = rand::thread_rng();
                let mut col = Vec3::new(0.0, 0.0, 0.0);
                for _ in 0..image_params.samples {
                    let u = (f64::from(i) + rng.gen::<f64>()) / f64::from(image_params.width);
                    let v = (f64::from(j) + rng.gen::<f64>()) / f64::from(image_params.height);
                    let ray = camera.get_ray(u, v);
                    col += color(&ray, world);
                }
                col = col / f64::from(image_params.samples);
                Vec3::new(col.x.sqrt(), col.y.sqrt(), col.z.sqrt())
            })
        })
        .collect();

    for color in colors {
        println!("{}", color.as_color_string())
    }
}
