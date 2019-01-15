use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, PartialEq, Copy, Clone)]
struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    fn r(&self) -> f64 {
        self.x
    }
    fn g(&self) -> f64 {
        self.y
    }
    fn b(&self) -> f64 {
        self.z
    }

    fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x, y, z }
    }
    fn squared_length(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    fn length(&self) -> f64 {
        self.squared_length().sqrt()
    }
    fn dot(&self, other: &Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    fn cross(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
    fn unit(&self) -> Vec3 {
        *self / self.length()
    }
}

impl From<(f64, f64, f64)> for Vec3 {
    fn from(tuple: (f64, f64, f64)) -> Vec3 {
        Vec3 {
            x: tuple.0,
            y: tuple.1,
            z: tuple.2,
        }
    }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f64) -> Vec3 {
        Vec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        rhs * self
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f64) -> Vec3 {
        Vec3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Div<Vec3> for f64 {
    type Output = Vec3;

    fn div(self, rhs: Vec3) -> Vec3 {
        rhs / self
    }
}

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
    fn hit_by_ray(&self, ray: &Ray) -> f64 {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * oc.dot(&ray.direction);
        let c = oc.dot(&oc) - self.radius * self.radius;
        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return -1.0;
        }
        (-b - discriminant.sqrt()) / (2.0 * a)
    }
}

fn color(ray: Ray) -> Vec3 {
    let sphere = Sphere::new(Vec3::new(0.0, 0.0, -1.0), 0.5);
    let t = sphere.hit_by_ray(&ray);
    if t > 0.0 {
        let normal = (ray.point_at_parameter(t) - Vec3::new(0.0, 0.0, -1.0)).unit();
        return 0.5 * Vec3::new(normal.x + 1.0, normal.y + 1.0, normal.z + 1.0);
    }
    let unit_direction = ray.direction.unit();
    let t = (unit_direction.y + 1.0) * 0.5;
    (1.0 - t) * Vec3::new(1.0, 1.0, 1.0) + t * Vec3::new(0.5, 0.7, 1.0)
}

fn main() {
    let nx = 200;
    let ny = 100;

    println!("P3");
    println!("{} {}", nx, ny);
    println!("255");

    let lower_left_corner = Vec3::new(-2.0, -1.0, -1.0);
    let horizontal = Vec3::new(4.0, 0.0, 0.0);
    let vertical = Vec3::new(0.0, 2.0, 0.0);
    let origin = Vec3::new(0.0, 0.0, 0.0);

    for j in (0..ny).rev() {
        for i in 0..nx {
            let u = f64::from(i) / f64::from(nx);
            let v = f64::from(j) / f64::from(ny);
            let ray = Ray::new(origin, lower_left_corner + u * horizontal + v * vertical);

            let col = color(ray);
            let ir = (255.99 * col.r()) as u8;
            let ig = (255.99 * col.g()) as u8;
            let ib = (255.99 * col.b()) as u8;
            println!("{} {} {}", ir, ig, ib);
        }
    }
}
