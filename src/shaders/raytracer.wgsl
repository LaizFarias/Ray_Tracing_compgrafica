const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f {
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray {
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera {
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn environment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f {
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}


fn check_ray_collision(r: ray, max: f32) -> hit_record {
  var spheresCount = i32(uniforms[19]);
  var quadsCount = i32(uniforms[20]);
  var boxesCount = i32(uniforms[21]);
  var meshCount = i32(uniforms[27]);

  var record = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);

  for (var i = 0; i < spheresCount; i++) {
    var record_sphere = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    var sphere_i = spheresb[i];
    var sphere_center = vec3f(sphere_i.transform[0], sphere_i.transform[1], sphere_i.transform[2]);
    var sphere_radius = f32(sphere_i.transform[3]);
    hit_sphere(sphere_center, sphere_radius, r, &record_sphere, max);

    if (record_sphere.hit_anything && record_sphere.t < record.t) {
      record = record_sphere;
      record.object_color = sphere_i.color;
      record.object_material = sphere_i.material;
    }
  }

  return record;
}

fn lambertian(normal: vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour {
  var rnd = rng_next_float(rng_state);

  if (rnd >= absorption) {
    var dir_result = normalize(normal + random_sphere);
    return material_behaviour(true, dir_result);
  }

  return material_behaviour(false, vec3f(0.0));
}

fn metal(normal: vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour {
  var reflected_vec = reflect(direction, normal);
  var final_dir = normalize(reflected_vec + random_sphere * fuzz);
  return material_behaviour(true, final_dir);
}

fn dielectric(normal: vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour { 
  var ratio = select(refraction_index, 1.0 / refraction_index, frontface);
  var unit_vec = normalize(r_direction);

  var cos_angle = min(dot(-unit_vec, normal), 1.0);
  var sin_angle = sqrt(1.0 - cos_angle * cos_angle);

  var cannot_refract = ratio * sin_angle > 1.0;
  var result_dir = vec3f(0.0);

  var base_reflect = pow((1.0 - refraction_index) / (1.0 + refraction_index), 2.0);
  var reflection_factor = base_reflect + (1.0 - base_reflect) * pow(1.0 - cos_angle, 5.0);
  var rand_factor = rng_next_float(rng_state);

  if (!(cannot_refract || rand_factor < reflection_factor)) {
    var perp = ratio * (unit_vec + cos_angle * normal);
    var parallel = -sqrt(abs(1.0 - dot(perp, perp))) * normal;
    result_dir = normalize(perp + parallel);
  } else {
    result_dir = reflect(unit_vec, normal);
  }

  return material_behaviour(true, result_dir);
}

fn emmisive(color: vec3f, light: f32) -> material_behaviour {
  return material_behaviour(false, vec3f(0.0));
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f {
  var limit = i32(uniforms[2]);
  var accumulated_light = vec3f(0.0);
  var attenuation = vec3f(1.0);
  var ray_temp = r;

  var bg_top = int_to_rgb(i32(uniforms[11]));
  var bg_bottom = int_to_rgb(i32(uniforms[12]));

  for (var bounce = 0; bounce < limit; bounce = bounce + 1) {
    var collision = check_ray_collision(ray_temp, RAY_TMAX);

    if (!collision.hit_anything) {
      accumulated_light += environment_color(ray_temp.direction, bg_top, bg_bottom) * attenuation;
      return accumulated_light;
    }

    var obj_color = vec3(collision.object_color[0], collision.object_color[1], collision.object_color[2]);
    var obj_mat = collision.object_material;

    var smoothness = obj_mat[0];
    var absorb = obj_mat[1];
    var fuzz = obj_mat[1];
    var specular = obj_mat[2];
    var ref_index = obj_mat[2];
    var emission = obj_mat[3];

    var rand_vec = rng_next_vec3_in_unit_sphere(rng_state);

    if (smoothness < 0.0) {
      var diel_behaviour = dielectric(collision.normal, ray_temp.direction, ref_index, collision.frontface, rand_vec, fuzz, rng_state);
      ray_temp.direction = diel_behaviour.direction;
      attenuation *= obj_color;
    } else {
      var lambert = lambertian(collision.normal, absorb, rand_vec, rng_state);
      var metal_behaviour = metal(collision.normal, ray_temp.direction, fuzz, rand_vec);

      if (specular != 0.0) {
        var rand_spec = rng_next_float(rng_state);
        if (rand_spec <= specular) {
          ray_temp.direction = metal_behaviour.direction;
          attenuation *= smoothness * vec3(1.0) + (1.0 - smoothness) * obj_color;
        } else {
          attenuation *= obj_color;
          if (lambert.scatter) {
            ray_temp.direction = lambert.direction;
          } else {
            accumulated_light += environment_color(ray_temp.direction, bg_top, bg_bottom) * attenuation;
            return accumulated_light;
          }
        }
      } else {
        if (smoothness > 0.0) {
          ray_temp.direction = metal_behaviour.direction;
          attenuation *= smoothness * vec3(1.0) + (1.0 - smoothness) * obj_color;
        } else {
          attenuation *= obj_color;
          if (lambert.scatter) {
            ray_temp.direction = lambert.direction;
          } else {
            accumulated_light += environment_color(ray_temp.direction, bg_top, bg_bottom) * attenuation;
            return accumulated_light;
          }
        }
      }
    }

    if (emission > 0.0) {
      accumulated_light += emission * obj_color * attenuation;
    }

    ray_temp.origin = collision.p;
  }

  accumulated_light += environment_color(ray_temp.direction, bg_top, bg_bottom) * attenuation;
  return accumulated_light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // init_rng (random number generator) we pass the pixel position, resolution and frame
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(rez)), time);

    // Get uv
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(rez);

    // Camera
    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);

    // Get camera
    var cam = get_camera(lookfrom, lookat, vec3(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);
    var samples_per_pixel = i32(uniforms[4]);

    var color = vec3f(0.0);

    // Steps:
    // 1. Loop for each sample per pixel
    // 2. Get ray
    // 3. Call trace function
    // 4. Average the color

    for (var i = 0; i < samples_per_pixel; i++) {
      var ray = get_ray(cam, uv, &rng_state);
      var light = trace(ray, &rng_state);
      color += light;
    } 

    color = color/f32(samples_per_pixel);

    var color_out = vec4(linear_to_gamma(color), 1.0);
    var map_fb = mapfb(id.xy, rez);

    // 5. Accumulate the color
    var should_accumulate = uniforms[3];
    var accumulated = rtfb[map_fb] * should_accumulate + color_out;

    // Set the color to the framebuffer
    rtfb[map_fb] = accumulated;
    fb[map_fb] = accumulated/ rtfb[map_fb].w;

}