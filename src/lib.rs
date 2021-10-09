
pub mod core {

    pub mod any_map {
        use ::core::any::{TypeId, Any};
        use ::core::hash::{BuildHasherDefault, Hasher};
        use hashbrown::HashMap;

        #[derive(Default)]
        pub struct NoOpHasher {
            hash: u64,
        }

        impl Hasher for NoOpHasher {
            fn write_u64(&mut self, n: u64) {
                debug_assert_eq!(self.hash, 0, "write_u64 should only be called the first time!");
                self.hash = n;
            }

            fn write_u128(&mut self, n: u128) {
                debug_assert_eq!(self.hash, 0, "write_u128 should only be called the first time!");
                self.hash = n as u64;
            }

            fn write(&mut self, _bytes: &[u8]) {
                panic!("NoOpHasher can only work on already hashed u64 (or u128) values!");
            }

            #[inline]
            fn finish(&self) -> u64 {
                // return cached hash
                self.hash
            }
        }

        /// HashMap meant to be used with TypeId keys (already fully hashed u64).
        pub type AnyMap = HashMap<TypeId, Box<dyn Any>, BuildHasherDefault<NoOpHasher>>;
    }

    pub mod math {
        use std::f32::consts::PI;
        use std::ops::{
            Neg,
            Add, AddAssign,
            Sub, SubAssign,
            Mul, MulAssign,
            Div, DivAssign
        };
        
        /// [`Vec2`]
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub struct Vec2 {
            pub x: f32,
            pub y: f32,
        }
        
        /// `vec2(x, y)' is a shortand for `Vec2 { x, y }`
        #[inline]
        pub fn vec2(x: f32, y: f32) -> Vec2 {
            Vec2 { x, y }
        }
        
        impl Vec2 {
            pub const POS_X: Self = Self { x: 1.0, y: 0.0 };  // Mask X or East direction
            pub const POS_Y: Self = Self { x: 0.0, y: 1.0 };  // Mask Y or South direction
            pub const NEG_X: Self = Self { x:-1.0, y: 0.0 };
            pub const NEG_Y: Self = Self { x: 0.0, y:-1.0 };
        
            pub const INVERT_X: Self = Self { x:-1.0, y: 1.0 };
            pub const INVERT_Y: Self = Self { x: 1.0, y:-1.0 };
        
            pub const ZERO: Self = Self { x: 0.0, y: 0.0 };
            pub const ONE: Self = Self { x: 1.0, y: 1.0 };
        
            pub const MIN: Self = Self { x: f32::MIN, y: f32::MIN };
            pub const MAX: Self = Self { x: f32::MAX, y: f32::MAX };
        
            #[inline]
            pub fn new(x: f32, y: f32) -> Self {
                Self { x, y }
            }
        
            /// Create a vector with the same x and y value.
            #[inline]
            pub fn full(v: f32) -> Self {
                Self { x: v, y: v }
            }

            /// Set x and y to the same value.
            #[inline]
            pub fn fill(&mut self, v: f32) -> &mut Self {
                self.x = v;
                self.y = v;
                self
            }

            /// Set x and y to zero.
            #[inline]
            pub fn zero(&mut self) -> &mut Self {
                self.x = 0.0;
                self.y = 0.0;
                self
            }
        
            /// Create a vector from length and angle (in radians).
            #[inline]
            pub fn polar(length: f32, angle: f32) -> Self {
                let (sin, cos) = angle.sin_cos();
        
                Self {
                    x: length * cos,
                    y: length * sin,
                }
            }
        
            /// Vec2 with unit norm (robust to zero norm).
            #[inline]
            pub fn normalize_or_zero(&self) -> Self {
                let norm = self.length();
        
                if norm > 0.0 {
                    Self { x: self.x / norm, y: self.y / norm }
                } else {
                    Self::ZERO
                }
            }
        
            /// Vec2 with unit norm (panics on zero norm).
            #[inline]
            pub fn normalize(&self) -> Self {
                let norm = self.length();
                Self { x: self.x / norm, y: self.y / norm }
            }
        
            /// The l2 norm (euclidean distance) of the vector.
            #[inline]
            pub fn length(&self) -> f32 {
                self.y.hypot(self.x)
            }
        
            /// The l2 norm (euclidean distance) of the vector, squared.
            #[inline]
            pub fn length_sq(&self) -> f32 {
                self.x * self.x + self.y * self.y
            }
        
            /// Sum the x and y component.
            #[inline]
            pub fn sum(&self) -> f32 {
                self.x + self.y
            }

            /// Vec2 normal (perpendicular) to the vector.
            #[inline]
            pub fn normal(&self) -> Self {
                let norm = self.length();
        
                if norm > 0.0 {
                    Self { x: -self.y / norm, y: self.x / norm }
                } else {
                    Self::ZERO
                }
            }
        
            /// Measures the angle of the vector (in radians).
            #[inline]
            pub fn angle(&self) -> f32 {
                self.y.atan2(self.x)
            }
        
            /// Rotate the vector by a given angle (in radians) w.r.t. the origin (or use `(self - pivot).rotate(angle) + pivot` to rotate w.r.t. a pivot point).
            #[inline]
            pub fn rotate(&self, angle: f32) -> Self {
                let (sin, cos) = angle.sin_cos();
        
                Self {
                    x: self.x * cos - self.y * sin,
                    y: self.x * sin + self.y * cos,
                }
            }

            /// Rotate the vector by a given angle (cos, sin) w.r.t. the origin.
            #[inline]
            pub fn rotate_cos_sin(&self, cos: f32, sin: f32) -> Self {        
                Self {
                    x: self.x * cos - self.y * sin,
                    y: self.x * sin + self.y * cos,
                }
            }
        
            /// Vec2 with min x and min y among the two vectors.
            #[inline]
            pub fn min(&self, other: Self) -> Self {
                Self {
                    x: self.x.min(other.x),
                    y: self.y.min(other.y),
                }
            }
        
            /// Vec2 with max x and max y among the two vectors.
            #[inline]
            pub fn max(&self, other: Self) -> Self {
                Self {
                    x: self.x.max(other.x),
                    y: self.y.max(other.y),
                }
            }
        
            /// Vec2 with x and y clamped in min-max range.
            #[inline]
            pub fn clamp(&self, min: Self, max: Self) -> Self {
                Self {
                    x: self.x.clamp(min.x, max.x),
                    y: self.y.clamp(min.y, max.y),
                }
            }
        
            /// Vec2 with x.floor() and y.floor().
            #[inline]
            pub fn floor(&self) -> Self {
                Self {
                    x: self.x.floor(),
                    y: self.y.floor(),
                }
            }
        
            /// Vec2 with x.round() and y.round().
            #[inline]
            pub fn round(&self) -> Self {
                Self {
                    x: self.x.round(),
                    y: self.y.round(),
                }
            }
        
            /// Vec2 with x.ceil() and y.ceil().
            #[inline]
            pub fn ceil(&self) -> Self {
                Self {
                    x: self.x.ceil(),
                    y: self.y.ceil(),
                }
            }
        
            /// Vec2 with x.fract() and y.fract().
            #[inline]
            pub fn fract(&self) -> Self {
                Self {
                    x: self.x.fract(),
                    y: self.y.fract(),
                }
            }
        
            /// Vec2 with x.abs() and y.abs().
            #[inline]
            pub fn abs(&self) -> Self {
                Self {
                    x: self.x.abs(),
                    y: self.y.abs(),
                }
            }
        
            /// Vec2 with x.signum() and y.signum().
            #[inline]
            pub fn signum(&self) -> Self {
                Self {
                    x: self.x.signum(),
                    y: self.y.signum(),
                }
            }
        
            /// Vec2 with e^x and e^y.
            #[inline]
            pub fn exp(&self) -> Self {
                Self {
                    x: self.x.exp(),
                    y: self.y.exp(),
                }
            }
        
            /// Vec2 with x^f.x and y^f.y.
            #[inline]
            pub fn powf(&self, f: Self) -> Self {
                Self {
                    x: self.x.powf(f.x),
                    y: self.y.powf(f.y),
                }
            }
        
        }
        
        // display
        impl std::fmt::Display for Vec2 {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "Vec2({}, {})", self.x, self.y)
            }
        }
        
        // back and forth from [f32; 2]
        impl From<[f32; 2]> for Vec2 {
            #[inline]
            fn from(v: [f32; 2]) -> Self {
                Self { x: v[0], y: v[1] }
            }
        }
        
        impl From<&[f32; 2]> for Vec2 {
            #[inline]
            fn from(v: &[f32; 2]) -> Self {
                Self { x: v[0], y: v[1] }
            }
        }
        
        impl From<Vec2> for [f32; 2] {
            #[inline]
            fn from(v: Vec2) -> Self {
                [v.x, v.y]
            }
        }
        
        impl From<&Vec2> for [f32; 2] {
            #[inline]
            fn from(v: &Vec2) -> Self {
                [v.x, v.y]
            }
        }
        
        // back and forth from (f32, f32)
        impl From<(f32, f32)> for Vec2 {
            #[inline]
            fn from(v: (f32, f32)) -> Self {
                Self { x: v.0, y: v.1 }
            }
        }
        
        impl From<&(f32, f32)> for Vec2 {
            #[inline]
            fn from(v: &(f32, f32)) -> Self {
                Self { x: v.0, y: v.1 }
            }
        }
        
        impl From<Vec2> for (f32, f32) {
            #[inline]
            fn from(v: Vec2) -> Self {
                (v.x, v.y)
            }
        }
        
        impl From<&Vec2> for (f32, f32) {
            #[inline]
            fn from(v: &Vec2) -> Self {
                (v.x, v.y)
            }
        }
        
        // ops
        impl Eq for Vec2 {}

        // neg
        impl Neg for Vec2 {
            type Output = Self;
        
            #[inline]
            fn neg(self) -> Self::Output {
                Self::Output {
                    x: -self.x,
                    y: -self.y,
                }
            }
        }
        
        // add
        impl Add for Vec2 {
            type Output = Self;
        
            #[inline]
            fn add(self, rhs: Self) -> Self::Output {
                Self::Output {
                    x: self.x + rhs.x,
                    y: self.y + rhs.y,
                }
            }
        }
        
        impl AddAssign for Vec2 {
            #[inline]
            fn add_assign(&mut self, rhs: Self) {
                self.x += rhs.x;
                self.y += rhs.y;
            }
        }
        
        // sub
        impl Sub for Vec2 {
            type Output = Self;
        
            #[inline]
            fn sub(self, rhs: Self) -> Self::Output {
                Self::Output {
                    x: self.x - rhs.x,
                    y: self.y - rhs.y,
                }
            }
        }
        
        impl SubAssign for Vec2 {
            #[inline]
            fn sub_assign(&mut self, rhs: Self) {
                self.x -= rhs.x;
                self.y -= rhs.y;
            }
        }
        
        // mul
        impl Mul<Vec2> for Vec2 {
            type Output = Self;
        
            #[inline]
            fn mul(self, rhs: Self) -> Self::Output {
                Self::Output {
                    x: self.x * rhs.x,
                    y: self.y * rhs.y,
                }
            }
        }
        
        impl Mul<f32> for Vec2 {
            type Output = Self;
        
            #[inline]
            fn mul(self, rhs: f32) -> Self::Output {
                Self::Output {
                    x: self.x * rhs,
                    y: self.y * rhs,
                }
            }
        }
        
        impl MulAssign<Vec2> for Vec2 {
            #[inline]
            fn mul_assign(&mut self, rhs: Self) {
                self.x *= rhs.x;
                self.y *= rhs.y;
            }
        }
        
        impl MulAssign<f32> for Vec2 {
            #[inline]
            fn mul_assign(&mut self, rhs: f32) {
                self.x *= rhs;
                self.y *= rhs;
            }
        }
        
        // rmul (i.e. f32 * Vec2)
        impl Mul<Vec2> for f32 {
            type Output = Vec2;
        
            #[inline]
            fn mul(self, rhs: Vec2) -> Self::Output {
                Self::Output {
                    x: self * rhs.x,
                    y: self * rhs.y,
                }
            }
        }
        
        // div
        impl Div<Vec2> for Vec2 {
            type Output = Self;
        
            #[inline]
            fn div(self, rhs: Self) -> Self::Output {
                Self::Output {
                    x: self.x / rhs.x,
                    y: self.y / rhs.y,
                }
            }
        }
        
        impl Div<f32> for Vec2 {
            type Output = Self;
        
            #[inline]
            fn div(self, rhs: f32) -> Self::Output {
                Vec2 {
                    x: self.x / rhs,
                    y: self.y / rhs,
                }
            }
        }
        
        impl DivAssign<Vec2> for Vec2 {
            #[inline]
            fn div_assign(&mut self, rhs: Self) {
                self.x /= rhs.x;
                self.y /= rhs.y;
            }
        }
        
        impl DivAssign<f32> for Vec2 {
            #[inline]
            fn div_assign(&mut self, rhs: f32) {
                self.x /= rhs;
                self.y /= rhs;
            }
        }
        
        // rdiv (i.e. f32 / Vec2)
        impl Div<Vec2> for f32 {
            type Output = Vec2;
        
            #[inline]
            fn div(self, rhs: Vec2) -> Self::Output {
                Self::Output {
                    x: self / rhs.x,
                    y: self / rhs.y,
                }
            }
        }

        /// [`Rect`]
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub struct Rect {
            pub min: Vec2,
            pub max: Vec2,
        }

        impl Rect {

            pub const ZERO: Self = Self { min: Vec2::ZERO, max: Vec2::ZERO };
            pub const ONE: Self = Self { min: Vec2::ONE, max: Vec2::ONE };

            pub const MIN: Self = Self { min: Vec2::MAX, max: Vec2::MIN };
            pub const MAX: Self = Self { min: Vec2::MIN, max: Vec2::MAX };

            #[inline]
            pub fn from_ltrb(l: f32, t: f32, r: f32, b: f32) -> Self {
                Self { min: vec2(l, t), max: vec2(r, b) }
            }

            #[inline]
            pub fn from_lt_rb(lt: Vec2, rb: Vec2) -> Self {
                Self { min: lt, max: rb }
            }

            #[inline]
            pub fn from_xywh(x: f32, y: f32, w: f32, h: f32) -> Self {
                Self { min: vec2(x, y), max: vec2(x + w, y + h) }
            }

            #[inline]
            pub fn from_xy_wh(xy: Vec2, wh: Vec2) -> Self {
                Self { min: xy, max: xy + wh }
            }

            /// Translate the rect by the given [`Vec2`] translation delta.
            #[inline]
            pub fn translate(self, delta: Vec2) -> Self {
                Self {
                    min: self.min + delta,
                    max: self.max + delta,
                }
            }

            /// Scale the rect by the given [`Vec2`] scale factor (using `min` as the pivot).
            #[inline]
            pub fn scale(self, factor: Vec2) -> Self {
                Self {
                    min: self.min,
                    max: self.min + (self.max - self.min) * factor
                }
            }

            /// Rotate the rect by the given angle (in radians) w.r.t. a given `abs_pivot` returning the rotated rect vertices in `lt`, `rt`, `lb`, `rb` order.
            #[inline]
            pub fn rotate(self, angle: f32, abs_pivot: Vec2) -> [Vec2; 4] {
                let (sin, cos) = angle.sin_cos();

                self.rotate_cos_sin(cos, sin, abs_pivot)
            }

            /// Rotate the rect by the given angle cos and sin w.r.t. a given `abs_pivot` returning the rotated rect vertices in `lt`, `rt`, `lb`, `rb` order.
            #[inline]
            pub fn rotate_cos_sin(self, cos: f32, sin: f32, abs_pivot: Vec2) -> [Vec2; 4] {
                [
                    abs_pivot + (self.lt() - abs_pivot).rotate_cos_sin(cos, sin),
                    abs_pivot + (self.rt() - abs_pivot).rotate_cos_sin(cos, sin),
                    abs_pivot + (self.lb() - abs_pivot).rotate_cos_sin(cos, sin),
                    abs_pivot + (self.rb() - abs_pivot).rotate_cos_sin(cos, sin),
                ]
            }

            /// Return the rectangular intersection of the two rects if they are colliding.
            #[inline]
            pub fn intersects(self, rect: Self) -> Option<Self> {
                let intersection = Self {
                    min: self.min.max(rect.min),
                    max: self.max.min(rect.max),
                };

                if intersection.is_neg() {
                    Some(intersection)
                } else {
                    None
                }
            }

            /// Return whether the point is inside the rect.
            #[inline]
            pub fn contains(&self, p: Vec2) -> bool {
                self.min.x <= p.x && p.x <= self.max.x && self.min.y <= p.y && p.y <= self.max.y
            }

            /// Return the left of the rect.
            #[inline]
            pub fn l(&self) -> f32 {
                self.min.x
            }

            /// Return the top of the rect.
            #[inline]
            pub fn t(&self) -> f32 {
                self.min.y
            }

            /// Return the right of the rect.
            #[inline]
            pub fn r(&self) -> f32 {
                self.max.x
            }

            /// Return the bottom of the rect.
            #[inline]
            pub fn b(&self) -> f32 {
                self.max.y
            }

            /// Return the topleft as a Vec2.
            #[inline]
            pub fn lt(&self) -> Vec2 {
                vec2(self.min.x, self.min.y)
            }

            /// Return the topright as a Vec2.
            #[inline]
            pub fn rt(&self) -> Vec2 {
                vec2(self.max.x, self.min.y) 
            }

            /// Return the bottomleft as a Vec2.
            #[inline]
            pub fn lb(&self) -> Vec2 {
                vec2(self.min.x, self.max.y)
            }

            /// Return the bottomright as a Vec2.
            #[inline]
            pub fn rb(&self) -> Vec2 {
                vec2(self.max.x, self.max.y)
            }

            /// Return [l, t, r, b].
            #[inline]
            pub fn ltrb(&self) -> [f32; 4] {
                [self.min.x, self.min.y, self.max.x, self.max.y]
            }

            /// Return (lt, rb) as a (Vec2, Vec2).
            #[inline]
            pub fn lt_rb(&self) -> (Vec2, Vec2) {
                (self.min, self.max)
            }

            /// Return the left of the rect.
            #[inline]
            pub fn x(&self) -> f32 {
                self.min.x
            }

            /// Return the top of the rect.
            #[inline]
            pub fn y(&self) -> f32 {
                self.min.y
            }

            // Return the topleft as a Vec2.
            #[inline]
            pub fn xy(&self) -> Vec2 {
                vec2(self.min.x, self.min.y)
            }

            /// Return the width of the rect.
            #[inline]
            pub fn w(&self) -> f32 {
                self.max.x - self.min.x
            }

            /// Return the height of the rect.
            #[inline]
            pub fn h(&self) -> f32 {
                self.max.y - self.min.y
            }

            /// Return the size of the rect.
            #[inline]
            pub fn wh(&self) -> Vec2 {
                self.max - self.min
            }

            /// Return the width of the rect.
            #[inline]
            pub fn width(&self) -> f32 {
                self.max.x - self.min.x
            }

            /// Return the height of the rect.
            #[inline]
            pub fn height(&self) -> f32 {
                self.max.y - self.min.y
            }
                
            /// Return the size of the rect.
            #[inline]
            pub fn size(&self) -> Vec2 {
                self.max - self.min
            }

            /// Return [x, y, w, h].
            #[inline]
            pub fn xywh(&self) -> [f32; 4] {
                [self.min.x, self.min.y, self.max.x - self.min.x, self.max.y - self.min.y]
            }

            /// Return (xy, wh) as a (Vec2, Vec2).
            #[inline]
            pub fn xy_wh(&self) -> (Vec2, Vec2) {
                (self.min, self.max - self.min)
            }

            /// Return center
            #[inline]
            pub fn center(&self) -> Vec2 {
                (self.min + self.max) * 0.5
            }

            /// Return mid x
            #[inline]
            pub fn midx(&self) -> f32 {
                (self.min.x + self.max.x) * 0.5
            }

            /// Return mid y
            #[inline]
            pub fn midy(&self) -> f32 {
                (self.min.y + self.max.y) * 0.5
            }

            /// Return midleft
            #[inline]
            pub fn midl(&self) -> Vec2 {
                vec2(self.l(), self.midy())
            }

            /// Return midtop
            #[inline]
            pub fn midt(&self) -> Vec2 {
                vec2(self.midx(), self.t())
            }

            /// Return midright
            #[inline]
            pub fn midr(&self) -> Vec2 {
                vec2(self.r(), self.midy())
            }

            /// Return midbottom
            #[inline]
            pub fn midb(&self) -> Vec2 {
                vec2(self.midx(), self.b())
            }

            /// Return the area of the rect
            #[inline]
            pub fn area(&self) -> f32 {
                self.w() * self.h()
            }

            /// Return true if `rect.w() < 0.0 || rect.h() < 0.0`.
            #[inline]
            pub fn is_neg(&self) -> bool {
                self.max.x < self.min.x || self.max.y < self.min.y
            }

            /// Return true if `rect.w() > 0.0 && rect.h() > 0.0`.
            #[inline]
            pub fn is_pos(&self) -> bool {
                self.min.x < self.max.x && self.min.y < self.max.y
            }

            /// Return the rect in normalized device coordinates.
            #[inline]
            pub fn ndc(&self, viewport_size: Vec2) -> Rect {
                let half_viewport_size = viewport_size * 0.5;  // center
                let factor_to_ndc = vec2(2.0 / viewport_size.x, -2.0 / viewport_size.y);

                Rect::from_lt_rb(
                (self.min - half_viewport_size) * factor_to_ndc,
                (self.max - half_viewport_size) * factor_to_ndc,
                )
            }

        }

        // display
        impl std::fmt::Display for Rect {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "Rect(x={}, y={}, w={}, h={})", self.x(), self.y(), self.w(), self.h())
            }
        }

        /// [`Rot2`]
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub struct Rot2 {
            /// `rel_pivot` to be offset by object origin and scaled by object size (e.g. (0.5, 0.5) is object center). 
            pub rel_pivot: Vec2,
            /// `angle` in radians to rotate from (wrt pivot).
            pub angle: f32,
        }

        impl Rot2 {
            pub const DEG_TO_RAD: f32 = PI / 180.0;
            pub const RAD_TO_DEG: f32 = 180.0 / PI;

            #[inline]
            pub fn new(rel_pivot: Vec2, angle: f32) -> Self {
                Self { rel_pivot, angle }
            }

            #[inline]
            pub fn abs_pivot(&self, rect: Rect) -> Vec2 {
                // any object has a min bounding box (axis aligned) wrapping it.
                rect.lt() + self.rel_pivot * (rect.rb() - rect.lt())
            } 

            #[inline]
            pub fn from_center(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(0.5, 0.5),
                    angle,
                }
            }

            #[inline]
            pub fn from_lb(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(0.0, 0.0),
                    angle,
                }
            }

            #[inline]
            pub fn from_rb(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(1.0, 0.0),
                    angle,
                }
            }
            
            #[inline]
            pub fn from_lt(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(0.0, 1.0),
                    angle,
                }
            }

            #[inline]
            pub fn from_rt(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(1.0, 1.0),
                    angle,
                }
            }

            #[inline]
            pub fn from_midl(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(0.0, 0.5),
                    angle,
                }
            }

            #[inline]
            pub fn from_midt(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(0.5, 1.0),
                    angle,
                }
            }
            
            #[inline]
            pub fn from_midr(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(1.0, 0.5),
                    angle,
                }
            }

            #[inline]
            pub fn from_midb(angle: f32) -> Self {
                Self {
                    rel_pivot: vec2(0.5, 0.0),
                    angle,
                }
            }

        }

        impl Eq for Rot2 {}
    }

    pub mod color {
        use std::ops::{
            Neg,
            Add, AddAssign,
            Sub, SubAssign,
            Mul, MulAssign,
            Div, DivAssign
        };
        
        /// [`Color`]
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub struct Color {
            pub r: f32,
            pub g: f32,
            pub b: f32,
            pub a: f32,
        }
        
        impl Color {
        
            // masks
            pub const R: Color = Color { r: 1.0, g: 0.0, b: 0.0, a: 0.0 };
            pub const G: Color = Color { r: 0.0, g: 1.0, b: 0.0, a: 0.0 };
            pub const B: Color = Color { r: 0.0, g: 0.0, b: 1.0, a: 0.0 };
            pub const A: Color = Color { r: 0.0, g: 0.0, b: 0.0, a: 1.0 };
        
            // constants
            pub const TRANSPARENT: Color = Color { r: 0.0, g: 0.0, b: 0.0, a: 0.0 };
            pub const BLACK: Color = Color { r: 0.0, g: 0.0, b: 0.0, a: 1.0 };
            pub const WHITE: Color = Color { r: 1.0, g: 1.0, b: 1.0, a: 1.0 };
            pub const RED: Color = Color { r: 1.0, g: 0.0, b: 0.0, a: 1.0 };
            pub const GREEN: Color = Color { r: 0.0, g: 1.0, b: 0.0, a: 1.0 };
            pub const BLUE: Color = Color { r: 0.0, g: 0.0, b: 1.0, a: 1.0 };

            #[inline]
            pub const fn rgba(r: f32, g: f32, b: f32, a: f32) -> Self {
                Self { r, g, b, a }
            }
        
            #[inline]
            pub fn clamp(self, min: f32, max: f32) -> Self {
                Self::rgba(
                    self.r.clamp(min, max), 
                    self.g.clamp(min, max), 
                    self.b.clamp(min, max), 
                    self.a.clamp(min, max),
                )
            }

            #[inline]
            pub fn brighten(self, amt: f32) -> Self {
                let brightness = self.r + self.g + self.b;
                let target_brightness = (brightness + amt).clamp(0.0, 3.0);

                Self::rgba(
                    self.r * target_brightness / brightness, 
                    self.g * target_brightness / brightness, 
                    self.b * target_brightness / brightness,
                    self.a
                )
            }

        }
        
        // display (1 decimal precision)
        impl ::std::fmt::Display for Color {
            fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
                write!(f, "Color({:.1}, {:.1}, {:.1}, {:.1})", self.r, self.g, self.b, self.a)
            }
        }
        
        // back and forth from [f32; 4]
        impl From<[f32; 4]> for Color {
            #[inline]
            fn from(rgba: [f32; 4]) -> Self {
                Self::rgba(rgba[0], rgba[1], rgba[2], rgba[3])
            }
        }
        
        impl From<&[f32; 4]> for Color {
            #[inline]
            fn from(rgba: &[f32; 4]) -> Self {
                Self::rgba(rgba[0], rgba[1], rgba[2], rgba[3])
            }
        }
        
        impl From<Color> for [f32; 4] {
            #[inline]
            fn from(color: Color) -> Self {
                [color.r, color.g, color.b, color.a]
            }
        }
        
        impl From<&Color> for [f32; 4] {
            #[inline]
            fn from(color: &Color) -> Self {
                [color.r, color.g, color.b, color.a]
            }
        }
        
        // back and forth from [u8; 4]
        impl From<[u8; 4]> for Color {
            #[inline]
            fn from(rgba: [u8; 4]) -> Self {
                Self::rgba(rgba[0] as f32 / 255.0, rgba[1] as f32 / 255.0, rgba[2] as f32 / 255.0, rgba[3] as f32 / 255.0)
            }
        }
        
        impl From<&[u8; 4]> for Color {
            #[inline]
            fn from(rgba: &[u8; 4]) -> Self {
                Self::rgba(rgba[0] as f32 / 255.0, rgba[1] as f32 / 255.0, rgba[2] as f32 / 255.0, rgba[3] as f32 / 255.0)
            }
        }
        
        impl From<Color> for [u8; 4] {
            #[inline]
            fn from(color: Color) -> Self {
                [(color.r * 255.0) as u8, (color.g * 255.0) as u8, (color.b * 255.0) as u8, (color.a * 255.0) as u8]
            }
        }
        
        impl From<&Color> for [u8; 4] {
            #[inline]
            fn from(color: &Color) -> Self {
                [(color.r * 255.0) as u8, (color.g * 255.0) as u8, (color.b * 255.0) as u8, (color.a * 255.0) as u8]
            }
        }
        
        // back and forth from (f32, f32, f32, f32)
        impl From<(f32, f32, f32, f32)> for Color {
            #[inline]
            fn from(rgba: (f32, f32, f32, f32)) -> Self {
                Self::rgba(rgba.0, rgba.1, rgba.2, rgba.3)
            }
        }
        
        impl From<&(f32, f32, f32, f32)> for Color {
            #[inline]
            fn from(rgba: &(f32, f32, f32, f32)) -> Self {
                Self::rgba(rgba.0, rgba.1, rgba.2, rgba.3)
            }
        }
        
        impl From<Color> for (f32, f32, f32, f32) {
            #[inline]
            fn from(color: Color) -> Self {
                (color.r, color.g, color.b, color.a)
            }
        }
        
        impl From<&Color> for (f32, f32, f32, f32) {
            #[inline]
            fn from(color: &Color) -> Self {
                (color.r, color.g, color.b, color.a)
            }
        }
        
        // back and forth from (u8, u8, u8, u8)
        impl From<(u8, u8, u8, u8)> for Color {
            #[inline]
            fn from(rgba: (u8, u8, u8, u8)) -> Self {
                Self::rgba(rgba.0 as f32 / 255.0, rgba.1 as f32 / 255.0, rgba.2 as f32 / 255.0, rgba.3 as f32 / 255.0)
            }
        }
        
        impl From<&(u8, u8, u8, u8)> for Color {
            #[inline]
            fn from(rgba: &(u8, u8, u8, u8)) -> Self {
                Self::rgba(rgba.0 as f32 / 255.0, rgba.1 as f32 / 255.0, rgba.2 as f32 / 255.0, rgba.3 as f32 / 255.0)
            }
        }
        
        impl From<Color> for (u8, u8, u8, u8) {
            #[inline]
            fn from(color: Color) -> Self {
                ((color.r * 255.0) as u8, (color.g * 255.0) as u8, (color.b * 255.0) as u8, (color.a * 255.0) as u8)
            }
        }
        
        impl From<&Color> for (u8, u8, u8, u8) {
            #[inline]
            fn from(color: &Color) -> Self {
                ((color.r * 255.0) as u8, (color.g * 255.0) as u8, (color.b * 255.0) as u8, (color.a * 255.0) as u8)
            }
        }
        
        // ops
        impl Eq for Color {}
        
        // neg
        impl Neg for Color {
            type Output = Self;
        
            #[inline]
            fn neg(self) -> Self::Output {
                Self::Output {
                    r: -self.r,
                    g: -self.g,
                    b: -self.b,
                    a: -self.a,
                }
            }
        }
        
        // add
        impl Add<Color> for Color {
            type Output = Self;
        
            fn add(self, color: Color) -> Self::Output {
                Self::Output::rgba(self.r + color.r, self.g + color.g, self.b + color.b, self.a + color.a)
            }
        }
        
        impl Add<f32> for Color {
            type Output = Self;
        
            fn add(self, c: f32) -> Self::Output {
                Self::Output::rgba(self.r + c, self.g + c, self.b + c, self.a + c)
            }
        }
        
        impl AddAssign<Color> for Color {
            fn add_assign(&mut self, color: Color) {
                self.r += color.r;
                self.g += color.g;
                self.b += color.b;
                self.a += color.a;
            }
        }
        
        impl AddAssign<f32> for Color {
            fn add_assign(&mut self, c: f32) {
                self.r += c;
                self.g += c;
                self.b += c;
                self.a += c;
            }
        }
        
        // sub
        impl Sub<Color> for Color {
            type Output = Self;
        
            fn sub(self, color: Color) -> Self::Output {
                Self::Output::rgba(self.r - color.r, self.g - color.g, self.b - color.b, self.a - color.a)
            }
        }
        
        impl Sub<f32> for Color {
            type Output = Self;
        
            fn sub(self, c: f32) -> Self::Output {
                Self::Output::rgba(self.r - c, self.g - c, self.b - c, self.a - c)
            }
        }
        
        impl SubAssign<Color> for Color {
            fn sub_assign(&mut self, color: Color) {
                self.r -= color.r;
                self.g -= color.g;
                self.b -= color.b;
                self.a -= color.a;
            }
        }
        
        impl SubAssign<f32> for Color {
            fn sub_assign(&mut self, c: f32) {
                self.r -= c;
                self.g -= c;
                self.b -= c;
                self.a -= c;
            }
        }
        
        // mul
        impl Mul<Color> for Color {
            type Output = Self;
        
            fn mul(self, color: Color) -> Self::Output {
                Self::Output::rgba(self.r * color.r, self.g * color.g, self.b * color.b, self.a * color.a)
            }
        }
        
        impl Mul<f32> for Color {
            type Output = Self;
        
            fn mul(self, c: f32) -> Self::Output {
                Self::Output::rgba(self.r * c, self.g * c, self.b * c, self.a * c)
            }
        }
        
        impl MulAssign<Color> for Color {
            fn mul_assign(&mut self, color: Color) {
                self.r *= color.r;
                self.g *= color.g;
                self.b *= color.b;
                self.a *= color.a;
            }
        }
        
        impl MulAssign<f32> for Color {
            fn mul_assign(&mut self, c: f32) {
                self.r *= c;
                self.g *= c;
                self.b *= c;
                self.a *= c;
            }
        }
        
        // rmul (i.e. f32 * Color)
        impl Mul<Color> for f32 {
            type Output = Color;
        
            #[inline]
            fn mul(self, rhs: Color) -> Self::Output {
                Self::Output {
                    r: self * rhs.r,
                    g: self * rhs.g,
                    b: self * rhs.b,
                    a: self * rhs.a,
                }
            }
        }
        
        // div
        impl Div<Color> for Color {
            type Output = Self;
        
            fn div(self, color: Color) -> Self::Output {
                Self::Output::rgba(self.r / color.r, self.g / color.g, self.b / color.b, self.a / color.a)
            }
        }
        
        impl Div<f32> for Color {
            type Output = Self;
        
            fn div(self, c: f32) -> Self::Output {
                Self::Output::rgba(self.r / c, self.g / c, self.b / c, self.a / c)
            }
        }
        
        impl DivAssign<Color> for Color {
            fn div_assign(&mut self, color: Color) {
                self.r /= color.r;
                self.g /= color.g;
                self.b /= color.b;
                self.a /= color.a;
            }
        }
        
        impl DivAssign<f32> for Color {
            fn div_assign(&mut self, c: f32) {
                self.r /= c;
                self.g /= c;
                self.b /= c;
                self.a /= c;
            }
        }
        
        // rdiv (i.e. f32 / Color)
        impl Div<Color> for f32 {
            type Output = Color;
        
            #[inline]
            fn div(self, rhs: Color) -> Self::Output {
                Self::Output {
                    r: self / rhs.r,
                    g: self / rhs.g,
                    b: self / rhs.b,
                    a: self / rhs.a,
                }
            }
        }
        
    }

    pub mod draw {
        use crate::backend::Gpu;
        use wgpu::util::DeviceExt;
        use std::sync::Arc;
        use std::ops::Range;
        
        /// [`Vertex`] trait to describe vertex attributes and buffer layout.
        pub trait Vertex: Sized + Clone + bytemuck::Pod + bytemuck::Zeroable {
            const ATTRIBUTES: &'static [wgpu::VertexAttribute];
            const BUFFER_LAYOUT: wgpu::VertexBufferLayout<'static> = wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<Self>() as _,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: Self::ATTRIBUTES,
            };
        }
        
        /// [`UniformBlock`] trait to describe a shader's bind groups layout.
        pub trait UniformBlock: Sized + Clone + bytemuck::Pod + bytemuck::Zeroable {
            const UNIFORM_BIND_GROUP_LAYOUT_ENTRIES: &'static [wgpu::BindGroupLayoutEntry] = &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::all(),
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: std::num::NonZeroU64::new(std::mem::size_of::<Self>() as _),
                    },
                    count: None,
                },
            ];
        
            const TEXTURE_BIND_GROUP_LAYOUT_ENTRIES: &'static [wgpu::BindGroupLayoutEntry] = &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::all(),
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float{ filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::all(),
                    ty: wgpu::BindingType::Sampler {
                        comparison: false,
                        filtering: true,
                    },
                    count: None,
                },
            ];
        
            fn uniform_bind_group_entries<'a>(uniform_buffers: &'a [Arc<wgpu::Buffer>]) -> Vec<wgpu::BindGroupEntry<'a>> {
                let mut entries = Vec::with_capacity(Self::UNIFORM_BIND_GROUP_LAYOUT_ENTRIES.len());
        
                for layout_entry in Self::UNIFORM_BIND_GROUP_LAYOUT_ENTRIES.iter() {
                    entries.push(wgpu::BindGroupEntry {
                        binding: layout_entry.binding,
                        resource: wgpu::BindingResource::Buffer(
                            wgpu::BufferBinding {
                                buffer: &*uniform_buffers[layout_entry.binding as usize],
                                offset: 0,
                                size: None,
                            },
                        ),
                    });
                }
        
                entries
            }
        
            fn texture_bind_group_entries<'a>(texture_views_samplers: &'a [Arc<(wgpu::TextureView, wgpu::Sampler)>]) -> Vec<wgpu::BindGroupEntry<'a>> {
                let mut entries = Vec::with_capacity(Self::TEXTURE_BIND_GROUP_LAYOUT_ENTRIES.len());
        
                for (b, arc_view_sampler) in texture_views_samplers.iter().enumerate() {
                    entries.push(wgpu::BindGroupEntry {
                        binding: (2 * b) as u32,
                        resource: wgpu::BindingResource::TextureView(&arc_view_sampler.0),
                    });
        
                    entries.push(wgpu::BindGroupEntry {
                        binding: (2 * b + 1) as u32,
                        resource: wgpu::BindingResource::Sampler(&arc_view_sampler.1),
                    });
                }
        
                entries
            }
        
        }
        
        /// [`Texture`] wraps [`wgpu::Texture`] and the metadata needed to create a view and a sampler.
        #[derive(Debug)]
        pub struct Texture {
            pub raw: wgpu::Texture,
            pub extent: wgpu::Extent3d,
            pub format: wgpu::TextureFormat,
            pub address_mode: wgpu::AddressMode,
            pub filter_mode: wgpu::FilterMode,
        }
        
        impl Texture {

            pub fn from_bytes(
                gpu: Arc<Gpu>,
                data: &[u8],
                bpp: usize,
                width: usize,
                height: usize,
                format: wgpu::TextureFormat,
                address_mode: wgpu::AddressMode,
                filter_mode: wgpu::FilterMode,
            ) -> Self {
                let extent = wgpu::Extent3d { width: width as _, height: height as _, depth_or_array_layers: 1 };

                let wgpu_texture = gpu.device.create_texture(&wgpu::TextureDescriptor {
                    label: None,
                    size: extent,
                    mip_level_count: 1,
                    sample_count: 1,
                    dimension: wgpu::TextureDimension::D2,
                    format,
                    usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_DST | wgpu::TextureUsages::TEXTURE_BINDING,
                });

                gpu.queue.write_texture(
                    wgpu::ImageCopyTextureBase {
                        texture: &wgpu_texture,
                        mip_level: 0,
                        origin: wgpu::Origin3d::ZERO,
                        aspect: wgpu::TextureAspect::All,
                    },
                    data,
                    wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: std::num::NonZeroU32::new(bpp as u32 * width as u32),
                        rows_per_image: std::num::NonZeroU32::new(height as u32),
                    },
                    extent,
                );
            
                Self {
                    raw: wgpu_texture,
                    extent,
                    format,
                    address_mode,
                    filter_mode,
                }
            }

            #[inline]
            pub fn from_rgba8(gpu: Arc<Gpu>, data: &[u8], width: usize, height: usize) -> Self {
                debug_assert!(data.len() == (width * height * 4) as usize, "invalid dimensions!");
                
                Self::from_bytes(gpu, data, 4, width, height, wgpu::TextureFormat::Rgba8Unorm, wgpu::AddressMode::ClampToEdge, wgpu::FilterMode::Nearest)
            }

            #[inline]
            pub fn from_a8(gpu: Arc<Gpu>, data: &[u8], width: usize, height: usize) -> Self {
                debug_assert!(data.len() == (width * height) as usize, "invalid dimensions!");
            
                Self::from_bytes(gpu, data, 1, width, height, wgpu::TextureFormat::R8Unorm, wgpu::AddressMode::ClampToEdge, wgpu::FilterMode::Nearest)
            }

            pub fn from_rgba32(gpu: Arc<Gpu>, data: &[u32], width: usize, height: usize) -> Self {
                debug_assert!(data.len() == (width * height) as usize, "invalid dimensions!");

                let mut bytes = Vec::with_capacity(data.len() * 4);

                for rgba32 in data.iter().cloned() {
                    bytes.push((rgba32 & 0xFF) as u8);
                    bytes.push((rgba32 & 0x00FF) as u8);
                    bytes.push((rgba32 & 0x0000FF) as u8);
                    bytes.push((rgba32 & 0x000000FF) as u8);
                }
            
                Self::from_rgba8(gpu, &bytes, width, height)
            }
        
            pub fn load(
                gpu: Arc<Gpu>,
                path: &str,
                flip_x: bool,
                flip_y: bool,
            ) -> Self {
                use image::GenericImageView;
                
                let img = image::open(&std::path::Path::new(path))
                    .expect(&format!("unable to load texture from path {}", path));

                let width = img.width() as usize;
                let height = img.height() as usize;
                let data = img.into_bytes();

                let bpp = data.len() / (width as usize * height as usize);
                let stride = (width as usize * bpp) as usize;

                let data_vec = {
                    if flip_x && flip_y {
                        data.iter().rev().cloned().collect::<_>()
                    } else {
                        if flip_x {
                            let mut v = Vec::with_capacity(data.len());
                            for chunk in data.chunks_exact(stride) { v.extend(chunk.iter().rev()); }
                            v
                        } else if flip_y {
                            let mut v = Vec::with_capacity(data.len());
                            for chunk in data.chunks_exact(stride).rev() { v.extend(chunk); }
                            v
                        } else {
                            data.to_vec()
                        }
                    }
                };

                let format = match bpp {
                    1 => wgpu::TextureFormat::R8Unorm,
                    4 => wgpu::TextureFormat::Rgba8Unorm,
                    _ => panic!("bpp texture format not supported!"),
                };

                Self::from_bytes(
                    gpu,
                    &data_vec,
                    bpp,
                    width,
                    height, 
                    format,
                    wgpu::AddressMode::ClampToEdge,
                    wgpu::FilterMode::Linear,
                )
            }
        
            #[inline]
            pub fn w(&self) -> f32 {
                self.extent.width as _
            }
        
            #[inline]
            pub fn h(&self) -> f32 {
                self.extent.height as _
            }
        
            #[inline]
            pub fn width(&self) -> usize {
                self.extent.width as _
            }
        
            #[inline]
            pub fn height(&self) -> usize {
                self.extent.height as _
            }
        
            #[inline]
            pub fn depth(&self) -> usize {
                self.extent.depth_or_array_layers as _
            }

            #[inline]
            pub fn bpp(&self) -> usize {
                self.format.describe().block_size as _
            }

            #[inline]
            pub fn stride(&self) -> usize {
                (self.extent.width as u32 * self.format.describe().block_size as u32) as _
            }
        
        }

        /// save a bitmap to disk (`path` must include the extension).
        // ideally method on Texture which knows its data, bpp, w, h (TODO: how to retrieve tex data from gpu memory? TODO: flipping).
        pub fn save_bitmap(path: &str, data: &[u8], width: u16, height: u16) {
            let bpp = data.len() / (width as usize * height as usize);

            let w = width as u32;
            let h = height as u32;
            
            match bpp {
                1 => {
                    let mut im = image::DynamicImage::new_luma8(w, h).to_luma8();
                    for (i, c) in data.chunks_exact(bpp).enumerate() {
                        let y = i as u32 / w;
                        let x = i as u32 % w;
                        im.put_pixel(x, y, image::Luma([c[0]]));
                    }
                    im.save(path).unwrap();
                },
                3 => {
                    let mut im = image::DynamicImage::new_rgb8(w, h).to_rgb8();
                    for (i, c) in data.chunks_exact(bpp).enumerate() {
                        let y = i as u32 / w;
                        let x = i as u32 % w;
                        im.put_pixel(x, y, image::Rgb([c[0], c[1], c[2]]));
                    }
                    im.save(path).unwrap();
                },
                4 => {
                    let mut im = image::DynamicImage::new_rgba8(w, h).to_rgba8();
                    for (i, c) in data.chunks_exact(bpp).enumerate() {
                        let y = i as u32 / w;
                        let x = i as u32 % w;
                        im.put_pixel(x, y, image::Rgba([c[0], c[1], c[2], c[3]]));
                    }
                    im.save(path).unwrap();
                },
                _ => {
                    panic!("cannot save texture with {} channels: only A (1), RGB (3), RGBA (4) allowed.", bpp);
                }
            }
        }

        /// [`Shader`] describes the pipeline spec.
        #[derive(Debug)]
        pub struct Shader {
            pub gpu: Arc<Gpu>,
            pub module: wgpu::ShaderModule,
            pub pipeline: wgpu::RenderPipeline,
            pub pipeline_layout: wgpu::PipelineLayout,
            pub uniform_bind_group_layout: Option<wgpu::BindGroupLayout>, 
            pub texture_bind_group_layout: Option<wgpu::BindGroupLayout>,
            // we cannot afford to be generic over V, U so we keep these &'static here just to assert the V, U used as generic in the functions are compliant with the shader.
            pub vertex_buffer_layout: &'static wgpu::VertexBufferLayout<'static>,
            pub uniform_bind_group_layout_entries: &'static [wgpu::BindGroupLayoutEntry],
            pub texture_bind_group_layout_entries: &'static [wgpu::BindGroupLayoutEntry],
        }
        
        impl Shader {
        
            pub fn wgsl<V: Vertex, U: UniformBlock>(gpu: Arc<Gpu>, wgsl: &str) -> Self {
                let shader_module = gpu.device.create_shader_module(&wgpu::ShaderModuleDescriptor {
                    label: None,
                    source: wgpu::ShaderSource::Wgsl(wgsl.into()),
                });
        
                Self::new::<V, U>(gpu, shader_module)
            }
        
            pub fn new<V: Vertex, U: UniformBlock>(gpu: Arc<Gpu>, shader_module: wgpu::ShaderModule) -> Self {
                let uniform_bind_group_layout = {
                    if U::UNIFORM_BIND_GROUP_LAYOUT_ENTRIES.len() > 0 {
                        Some(gpu.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                            label: None,
                            entries: U::UNIFORM_BIND_GROUP_LAYOUT_ENTRIES,
                        }))
                    } else {
                        None
                    }
                };
        
                let texture_bind_group_layout = {
                    if U::TEXTURE_BIND_GROUP_LAYOUT_ENTRIES.len() > 0 {
                        Some(gpu.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                            label: None,
                            entries: U::TEXTURE_BIND_GROUP_LAYOUT_ENTRIES,
                        }))
                    } else {
                        None
                    }
                };
        
                let bind_group_layouts = {
                    if uniform_bind_group_layout.is_some() && texture_bind_group_layout.is_some() {
                        vec![uniform_bind_group_layout.as_ref().unwrap(), texture_bind_group_layout.as_ref().unwrap()]   
                    } else {
                        vec![uniform_bind_group_layout.as_ref().unwrap_or(texture_bind_group_layout.as_ref().unwrap())]
                    }
                };
        
                let pipeline_layout = gpu.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: None,
                    bind_group_layouts: &bind_group_layouts,
                    push_constant_ranges: &[],
                });
            
                let pipeline = gpu.device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                    label: None,
                    layout: Some(&pipeline_layout),
                    vertex: wgpu::VertexState {
                        module: &shader_module,
                        entry_point: "main",
                        buffers: &[V::BUFFER_LAYOUT],
                    },
                    fragment: Some(wgpu::FragmentState {
                        module: &shader_module,
                        entry_point: "main",
                        targets: &[wgpu::ColorTargetState {
                            format: wgpu::TextureFormat::Bgra8UnormSrgb,
                            blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                            write_mask: wgpu::ColorWrites::ALL,
                        }],
                    }),
                    primitive: wgpu::PrimitiveState {
                        topology: wgpu::PrimitiveTopology::TriangleList,
                        strip_index_format: None,
                        front_face: wgpu::FrontFace::Ccw,
                        cull_mode: None,
                        polygon_mode: wgpu::PolygonMode::Fill,
                        clamp_depth: false,
                        conservative: false,
                    },
                    depth_stencil: None,
                    multisample: wgpu::MultisampleState {
                        count: 1,
                        mask: !0,
                        alpha_to_coverage_enabled: false,
                    },
                });
        
                Self {
                    gpu,
                    module: shader_module,
                    pipeline,
                    pipeline_layout,
                    vertex_buffer_layout: &V::BUFFER_LAYOUT,
                    uniform_bind_group_layout_entries: U::UNIFORM_BIND_GROUP_LAYOUT_ENTRIES,
                    texture_bind_group_layout_entries: U::TEXTURE_BIND_GROUP_LAYOUT_ENTRIES,
                    uniform_bind_group_layout,
                    texture_bind_group_layout,
                }
            }
        
        }
        
        /// [`MakeShader`] trait.
        pub trait MakeShader {
        
            type Vertex: Vertex;
            type UniformBlock: UniformBlock;
        
            const SOURCE_WGSL: &'static str;
        
            fn make(gpu: Arc<Gpu>) -> Shader {
                Shader::wgsl::<Self::Vertex, Self::UniformBlock>(gpu, Self::SOURCE_WGSL)
            }
        
        }
        
        /// [`DrawOp`] wraps the data needed to issue a draw call.
        #[derive(Debug, Clone)]
        pub struct DrawOp {
            pub shader: Arc<Shader>,
            pub uniform_buffers: Option<Vec<Arc<wgpu::Buffer>>>,
            pub texture_views_samplers: Option<Vec<Arc<(wgpu::TextureView, wgpu::Sampler)>>>,
            pub uniform_bind_group: Option<Arc<wgpu::BindGroup>>,
            pub texture_bind_group: Option<Arc<wgpu::BindGroup>>,
            pub vertex_buffer: Option<Arc<wgpu::Buffer>>,
            pub index_buffer: Option<Arc<wgpu::Buffer>>,
            pub draw_range: Range<usize>,
        }
        
        impl DrawOp {
            // TODO: Ergonomics
            // - builder pattern
            // TODO: Performance
            // - update bindings in place (even by directly providing Arc<ResourceToShare>?)
            // - no heap allocation: SmallVec replacing Vec or fixed size array for unif buffers and textures
        
            pub fn new(shader: Arc<Shader>) -> Self {
                Self {
                    shader,
                    uniform_buffers: None,
                    texture_views_samplers: None,
                    uniform_bind_group: None,
                    texture_bind_group: None,
                    vertex_buffer: None,
                    index_buffer: None,
                    draw_range: 0 .. 0,
                }
            }
        
            pub fn bind_uniforms<U: UniformBlock>(&mut self, uniforms: Option<&[U]>, textures: Option<&[Arc<Texture>]>) {
                debug_assert_eq!(U::UNIFORM_BIND_GROUP_LAYOUT_ENTRIES, self.shader.uniform_bind_group_layout_entries, "UniformBlock 'uniform bind group layout entries' not compliant with Shader spec!");
                debug_assert_eq!(U::TEXTURE_BIND_GROUP_LAYOUT_ENTRIES, self.shader.texture_bind_group_layout_entries, "UniformBlock 'texture bind group layout entries' not compliant with Shader spec!");
        
                let gpu = &self.shader.gpu;

                if let Some(uniforms_) = uniforms {
                    // TODO: update existing if already initialized!
                    self.uniform_buffers = Some(uniforms_.iter().cloned().map(
                        |u| Arc::new(gpu.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                            label: None,
                            contents: bytemuck::cast_slice(&[u]),
                            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
                        }))
                    ).collect::<Vec<_>>());
        
                    self.uniform_bind_group = Some(Arc::new(gpu.device.create_bind_group(&wgpu::BindGroupDescriptor {
                        label: None,
                        layout: self.shader.uniform_bind_group_layout.as_ref().unwrap(),
                        entries: &U::uniform_bind_group_entries(&self.uniform_buffers.as_ref().unwrap()[..]),
                    })));
                }
        
                if let Some(textures_) = textures {
                    // TODO: allow sampler settings as additional parameters, consider not recreating the sampler every time.
                    // More in general, more efficient resource sharing (e.g. reusing the entire texture bind group).
                    self.texture_views_samplers = Some(textures_.iter().enumerate().map(
                        |(i, texture)| {
                            let view_dim = match U::TEXTURE_BIND_GROUP_LAYOUT_ENTRIES[2 * i].ty {
                                wgpu::BindingType::Texture { view_dimension, .. } => view_dimension,
                                _ => panic!("invalid texture bind group layout! (layout entry at {} should contain a texture view!)", 2 * i),
                            };
        
                            let view = texture.raw.create_view(&wgpu::TextureViewDescriptor {
                                label: None,
                                format: Some(texture.format),
                                dimension: Some(view_dim),
                                ..Default::default()
                            });
        
                            let sampler = gpu.device.create_sampler(&wgpu::SamplerDescriptor {
                                label: None,
                                address_mode_u: texture.address_mode,
                                address_mode_v: texture.address_mode,
                                address_mode_w: texture.address_mode,
                                mag_filter: texture.filter_mode,
                                min_filter: texture.filter_mode,
                                mipmap_filter: texture.filter_mode,
                                ..Default::default()
                            });
        
                            Arc::new((view, sampler))
                        }
                    ).collect::<Vec<_>>());
        
                    self.texture_bind_group = Some(Arc::new(gpu.device.create_bind_group(&wgpu::BindGroupDescriptor {
                        label: None,
                        layout: self.shader.texture_bind_group_layout.as_ref().unwrap(),
                        entries: &U::texture_bind_group_entries(&self.texture_views_samplers.as_ref().unwrap()[..]),
                    })));
                }
            }
        
            pub fn bind_vertices<V: Vertex>(&mut self, vertices: &[V]) {
                debug_assert_eq!(&V::BUFFER_LAYOUT, self.shader.vertex_buffer_layout, "Vertex 'buffer layout' not compliant with Shader spec!");
        
                let gpu = &self.shader.gpu;

                // TODO: update existing (from offset 0 or paramter) if already initialized!
                self.vertex_buffer = Some(Arc::new(gpu.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: None,
                    contents: bytemuck::cast_slice(vertices),
                    usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                })));
            }
        
            pub fn bind_indices(&mut self, indices: &[u16]) {
                let gpu = &self.shader.gpu;

                self.index_buffer = Some(Arc::new(gpu.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: None,
                    contents: bytemuck::cast_slice(indices),
                    usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
                })));
            }
        
        }

        /// [`prefabs`] ready-made [`Vertex`], [`UniformBlock`] and [`MakeShader`] recipes.
        pub mod prefabs {
            use super::{Vertex, UniformBlock, MakeShader};
        
            #[repr(C)]
            #[derive(Debug, Clone, Copy)]
            pub struct Vertex4f(pub [f32; 4]);
            unsafe impl bytemuck::Pod for Vertex4f {}
            unsafe impl bytemuck::Zeroable for Vertex4f {}
            impl Vertex for Vertex4f {
                const ATTRIBUTES: &'static [wgpu::VertexAttribute] = &wgpu::vertex_attr_array![
                    0 => Float32x4,
                ];
            }
            
            #[repr(C)]
            #[derive(Debug, Clone, Copy)]
            pub struct Vertex4fx2(pub [[f32; 4]; 2]);
            unsafe impl bytemuck::Pod for Vertex4fx2 {}
            unsafe impl bytemuck::Zeroable for Vertex4fx2 {}
            impl Vertex for Vertex4fx2 {
                const ATTRIBUTES: &'static [wgpu::VertexAttribute] = &wgpu::vertex_attr_array![
                    0 => Float32x4,
                    1 => Float32x4,
                ];
            }
        
            #[repr(C)]
            #[derive(Debug, Clone, Copy)]
            pub struct Sampler2d;
            unsafe impl bytemuck::Pod for Sampler2d {}
            unsafe impl bytemuck::Zeroable for Sampler2d {}
            impl UniformBlock for Sampler2d {
                const UNIFORM_BIND_GROUP_LAYOUT_ENTRIES: &'static [wgpu::BindGroupLayoutEntry] = &[];
            }

            #[derive(Debug)]
            pub struct MakeShader2d;

            impl MakeShader for MakeShader2d {

                type Vertex = Vertex4fx2;
                type UniformBlock = Sampler2d;

                const SOURCE_WGSL: &'static str = r#"
                struct VertexInput {
                    [[location(0)]] pos: vec4<f32>;
                    [[location(1)]] color: vec4<f32>;
                };
                
                struct VertexOutput {
                    [[builtin(position)]] clip_position: vec4<f32>;
                    [[location(0)]] uv: vec2<f32>;
                    [[location(1)]] color: vec4<f32>;
                };
                
                [[group(0), binding(0)]]
                var texture0: texture_2d<f32>;
                [[group(0), binding(1)]]
                var sampler: sampler;
                
                [[stage(vertex)]]
                fn main(in: VertexInput) -> VertexOutput {
                    var out: VertexOutput;
                    out.uv = in.pos.zw;
                    out.color = in.color;
                    out.clip_position = vec4<f32>(in.pos.xy, 0.0, 1.0);
                    return out;
                }
                
                [[stage(fragment)]]
                fn main(in: VertexOutput) -> [[location(0)]] vec4<f32> {
                    return in.color * textureSample(texture0, sampler, in.uv);
                }
                "#;

            }

            #[derive(Debug)]
            pub struct MakeShaderText2d;

            impl MakeShader for MakeShaderText2d {

                type Vertex = Vertex4fx2;
                type UniformBlock = Sampler2d;

                const SOURCE_WGSL: &'static str = r#"
                struct VertexInput {
                    [[location(0)]] pos: vec4<f32>;
                    [[location(1)]] color: vec4<f32>;
                };
                
                struct VertexOutput {
                    [[builtin(position)]] clip_position: vec4<f32>;
                    [[location(0)]] uv: vec2<f32>;
                    [[location(1)]] color: vec4<f32>;
                };
                
                [[group(0), binding(0)]]
                var texture0: texture_2d<f32>;
                [[group(0), binding(1)]]
                var sampler: sampler;
                
                [[stage(vertex)]]
                fn main(in: VertexInput) -> VertexOutput {
                    var out: VertexOutput;
                    out.uv = in.pos.zw;
                    out.color = in.color;
                    out.clip_position = vec4<f32>(in.pos.xy, 0.0, 1.0);
                    return out;
                }
                
                [[stage(fragment)]]
                fn main(in: VertexOutput) -> [[location(0)]] vec4<f32> {
                    return in.color * textureSample(texture0, sampler, in.uv).r;
                }
                "#;

            }

        }

    }

}


pub mod auxil {
    use crate::core::math::{Vec2, vec2};

    #[inline]
    pub fn ndc_pos(px_pos: Vec2, vp_size: Vec2) -> Vec2 {
        vec2(2.0, -2.0) * (px_pos - vp_size * 0.5) / vp_size
    }

    pub mod draw {
        use crate::core::math::{Vec2, vec2, Rect};
        use crate::core::color::Color;
        use crate::core::draw::{DrawOp, Vertex, UniformBlock, Texture, prefabs::{Vertex4fx2, Sampler2d}};
        use crate::auxil::draw::text::{CharDrawInfo, FontDrawInfo, TextStyle, TextDirection};
        use std::sync::Arc;

        pub mod text {
            use crate::core::math::{Vec2, Rect};
            use crate::core::color::Color;

            /// [`CharDrawInfo`] all info required to draw a rasterized char.
            #[derive(Debug, Clone, Copy)]
            pub struct CharDrawInfo {
                /// `px_region` is the texture rectangular region (in pixels) where this char is located.
                pub px_region: Rect,
                /// `px_advance` is the horizontal advance (in pixels) for this character.
                pub px_advance: f32,
            }

            /// [`FontDrawInfo`] all info required to draw text.
            pub trait FontDrawInfo {

                /// texture_size is needed to convert from pixels to normalized texture coords.
                fn texture_size(&self) -> Vec2;

                /// monospace fonts advance of a fixed amount (no matter the char pair).
                fn px_monospace_advance(&self) -> Option<f32> { None }
                /// proportional fonts advance of a variable amount read from a `kern_table[pre_char][cur_char]` (in addition to the width of the character).
                #[allow(unused_variables)]
                fn px_kern(&self, pre_char: char, cur_char: char) -> f32 { 0.0 }

                /// `CharDrawInfo` contains all the info required to draw a char (in particular its `px_region` in the texture and its `px_advance`). 
                fn get_char_draw_info_or_fallback(&mut self, ch: char) -> CharDrawInfo;

            }

            #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
            pub enum TextDirection {
                /// Left to Right.
                LTR,
                /// Right to Left.
                RTL,
            }

            impl Default for TextDirection {
                fn default() -> Self {
                    Self::LTR
                }
            }

            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct TextStyle {
                pub px_dest: Vec2,
                pub color: Color,
                pub direction: TextDirection,
                // scale, rotation, effects (bold, italic, underline, outline, glow)?
            }

            impl Default for TextStyle {
                fn default() -> Self {
                    Self {
                        px_dest: Vec2::ZERO,
                        color: Color::WHITE,
                        direction: Default::default(),
                    }
                }
            }

            // Monospace

            // (l2r: POS_X) ndc_rect l = dest.x + (fixed_advance.x - size.x) / 2
            // (r2l: NEG_X) ndc_rect l = dest.x - fixed_advance.x + (fixed_advance.x - size.x) / 2
            // (t2b: POS_Y) ndc_rect t = dest.y + (fixed_advance.y - size.y) / 2
            // (b2t: NEG_Y) ndc_rect t = dest.y - fixed_advance.y + (fixed_advance.y - size.y) / 2
            // (lt2rb: ONE) ndc_rect lt = dest + (fixed_advance - size) / 2
            // (rb2lt:-ONE) ndc_rect lt = dest - fixed_advance + (fixed_advance - size) / 2

            // (step: Vec2)
            // if step.x > 0.0 { step.x -= 1.0 };
            // if step.y > 0.0 { step.y -= 1.0 };
            // INCORRECT WIP!!!
            // let ndc_rect = Rect::from_xy_wh(
            //  dest + step * fixed_advance + step.signum() * (fixed_advance - size) / 2,
            //  size,
            // );
            
            // Proportional
            // 

        }

        /// [`Stream2d`] opinionistic attempt at batched primitive drawing.
        pub struct Stream2d {
            /// `draw_op` to be issued (contains the pipeline spec, the gpu buffers and textures).
            /// Most functions assume the draw_op's bound texture to have a white pixel at the very topleft (except when tex coords are provided directly).
            pub draw_op: DrawOp,
            /// `vertices` is a `Vec<Vertex4fx2>` whose internal layout is interpreted as `[[x, y, u, v], [r, g, b, a]]`.
            pub vertices: Vec<Vertex4fx2>,
            /// `indices` is a `Vec<u16>`.
            pub indices: Vec<u16>,
            /// `viewport` pixel size to handle conversion from/to normalized device coordinates. 
            pub viewport_size: Vec2,
            /// `is_indexed` drawing mode.
            pub is_indexed: bool,
            /// `font` to draw text.
            pub font: Option<Box<dyn FontDrawInfo>>,
        }

        impl Stream2d {

            pub fn new(draw_op: DrawOp) -> Self {
                debug_assert_eq!(draw_op.shader.uniform_bind_group_layout_entries, Sampler2d::UNIFORM_BIND_GROUP_LAYOUT_ENTRIES, "shader must support a Sampler2d UniformBlock");
                debug_assert_eq!(draw_op.shader.texture_bind_group_layout_entries, Sampler2d::TEXTURE_BIND_GROUP_LAYOUT_ENTRIES, "shader must support a Sampler2d UniformBlock");
                debug_assert_eq!(draw_op.shader.vertex_buffer_layout, &Vertex4fx2::BUFFER_LAYOUT, "shader must support a Vertex4fx2 Vertex!");

                Self {
                    draw_op,
                    vertices: Vec::new(),  // with_capacity(u16::MAX as usize) worth it? roughly 3 MB
                    indices: Vec::new(),
                    viewport_size: Vec2::ZERO,
                    is_indexed: false,
                    font: None,
                }
            }

            pub fn indexed(draw_op: DrawOp) -> Self {
                let mut pen = Self::new(draw_op);
                pen.is_indexed = true;
                pen
            } 

            /// clears the buffered vertices and indices. The gpu buffers, stored in the `DrawOp` are not cleared
            /// (they get updated only when binding, overwriting the relevant data portion, no smart caching atm!).
            pub fn clear(&mut self) {
                self.vertices.clear();
                self.indices.clear();
            }

            /// bind the current vertices and indices to the internal `DrawOp` updating the internal gpu buffers.
            pub fn bind_to_draw_op(&mut self) {
                self.draw_op.bind_vertices(&self.vertices);

                if self.indices.len() > 0 {
                    self.draw_op.bind_indices(&self.indices);
                    self.draw_op.draw_range = 0 .. self.indices.len();
                } else {
                    self.draw_op.index_buffer = None;
                    self.draw_op.draw_range = 0 .. self.vertices.len();
                }
            }

            /// bind the given font to this stream2d and the associated font atlas texture to the internal draw_op.
            pub fn bind_font<U: UniformBlock>(&mut self, font: Box<dyn FontDrawInfo>, uniforms: Option<&[U]>, textures: &[Arc<Texture>]) -> Option<Box<dyn FontDrawInfo>> {
                self.draw_op.bind_uniforms::<U>(uniforms, Some(textures));
                std::mem::replace(&mut self.font, Some(font))
            }

            /// draw a rect by providing a destination `px_rect` (in pixels) and a per-vertex `color` (e.g. fill).
            #[inline]
            pub fn fill_rect(&mut self, px_rect: Rect, color: Color) {
                self.fill_texture_nc(px_rect.ndc(self.viewport_size), Rect::ZERO, color);
            }

            /// draw a texture by providing a destination `ndc_rect` (in normalized device coordinates), a texture
            /// `ntc_rect` to sample from (in normalized texture coordinates) and a per-vertex `color` (e.g. blend rgba mul).
            #[inline]
            pub fn fill_texture_nc(&mut self, ndc_rect: Rect, ntc_rect: Rect, color: Color) {
                let rgba = color.into();

                if self.is_indexed {
                    let vtx = self.vertices.len() as u16;
                    let vertices = &mut self.vertices;
                    let indices = &mut self.indices;

                    vertices.push(Vertex4fx2([[ndc_rect.l(), ndc_rect.t(), ntc_rect.l(), ntc_rect.t()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.r(), ndc_rect.t(), ntc_rect.r(), ntc_rect.t()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.l(), ndc_rect.b(), ntc_rect.l(), ntc_rect.b()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.r(), ndc_rect.b(), ntc_rect.r(), ntc_rect.b()], rgba]));

                    indices.push(vtx); indices.push(vtx + 1); indices.push(vtx + 2);
                    indices.push(vtx + 1); indices.push(vtx + 2); indices.push(vtx + 3);
                } else {
                    let vertices = &mut self.vertices;
                    debug_assert!(self.indices.len() == 0, "cannot mix vertex-only with indexed calls!");

                    vertices.push(Vertex4fx2([[ndc_rect.l(), ndc_rect.t(), ntc_rect.l(), ntc_rect.t()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.r(), ndc_rect.t(), ntc_rect.r(), ntc_rect.t()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.l(), ndc_rect.b(), ntc_rect.l(), ntc_rect.b()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.r(), ndc_rect.t(), ntc_rect.r(), ntc_rect.t()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.l(), ndc_rect.b(), ntc_rect.l(), ntc_rect.b()], rgba]));
                    vertices.push(Vertex4fx2([[ndc_rect.r(), ndc_rect.b(), ntc_rect.r(), ntc_rect.b()], rgba]));
                }
            }

            /// draw a sequence of characters by providing the text `content_str` and `style`.
            #[inline]
            pub fn push_text(&mut self, content_str: &str, style: &TextStyle) {
                if let Some(mut font) = std::mem::take(&mut self.font) {
                    // assumes draw_op bound texture is a font atlas texture with metadata as described by this font!!!
                    let px_dest = style.px_dest;

                    let tex_size = font.texture_size();
                    let rec_tex_size = 1.0 / tex_size;
                    let vp_size = self.viewport_size;
                    let rec_vp_size_neg_y = Vec2::INVERT_Y / vp_size;

                    let mut ndc_dest = 2.0 * (px_dest - vp_size * 0.5) * rec_vp_size_neg_y;
    
                    if let Some(px_fixed_advance) = font.px_monospace_advance() {
                        let ndc_fixed_advance = px_fixed_advance * rec_vp_size_neg_y.x;
                        let ndc_fixed_advance_signed: f32 = match style.direction {
                            TextDirection::LTR => ndc_fixed_advance,
                            TextDirection::RTL => {
                                ndc_dest.x -= ndc_fixed_advance;
                                -ndc_fixed_advance
                            }
                        };
    
                        for ch in content_str.chars() {
                            let CharDrawInfo { px_region: ch_px_rect, px_advance: _ } = font.get_char_draw_info_or_fallback(ch);
                            let ch_ndc_size = ch_px_rect.size() * rec_vp_size_neg_y;
                            let ndc_rect = Rect::from_xy_wh(vec2(ndc_dest.x + (ndc_fixed_advance - ch_ndc_size.x) * 0.5, ndc_dest.y), ch_ndc_size);
                            let ntc_rect = Rect::from_xy_wh(ch_px_rect.xy() * rec_tex_size, ch_px_rect.wh() * rec_tex_size);
    
                            self.fill_texture_nc(ndc_rect, ntc_rect, style.color);
    
                            ndc_dest.x += ndc_fixed_advance_signed;
                        }
                    } else {
                        if let TextDirection::RTL = style.direction {
                            todo!("not implemented yet!");
                        }

                        let mut pre_ch = ' ';
    
                        for ch in content_str.chars() {
                            let CharDrawInfo { px_region: ch_px_rect, px_advance } = font.get_char_draw_info_or_fallback(ch);
                            let px_kern = font.px_kern(pre_ch, ch);
                            
                            ndc_dest.x += px_kern * rec_vp_size_neg_y.x;
    
                            if ch_px_rect.area() > 0.0 {  // if not space, tab or alike
                                let ch_ndc_size = ch_px_rect.size() * rec_vp_size_neg_y;
                                let ndc_rect = Rect::from_xy_wh(ndc_dest, ch_ndc_size);
                                let ntc_rect = Rect::from_xy_wh(ch_px_rect.xy() * rec_tex_size, ch_px_rect.wh() * rec_tex_size);
    
                                self.fill_texture_nc(ndc_rect, ntc_rect, style.color);
                            }
    
                            ndc_dest.x += px_advance * rec_vp_size_neg_y.x;
                            pre_ch = ch;
                        }
                    }
                
                    self.font = Some(font);
                } else {
                    eprintln!("Font not bound!");
                }
            }

        }

    }

    pub mod font {
        use crate::backend::Gpu;
        use crate::core::math::{Vec2, Rect};
        use crate::core::draw::Texture;
        use crate::auxil::draw::text::{CharDrawInfo, FontDrawInfo};
        use ab_glyph::{Font as abFont, ScaleFont as abScaleFont};
        use std::collections::HashMap;
        use std::sync::Arc;

        #[derive(Debug, Clone)]
        pub struct GlyphInfo {
            pub ab_id: ab_glyph::GlyphId,
            pub px_region: Rect,
            pub px_advance: f32,
        }

        impl From<&GlyphInfo> for CharDrawInfo {
            fn from(item: &GlyphInfo) -> Self {
                Self {
                    px_region: item.px_region.into(),
                    px_advance: item.px_advance as _,
                }
            }
        }

        impl From<GlyphInfo> for CharDrawInfo {
            fn from(item: GlyphInfo) -> Self {
                Self {
                    px_region: item.px_region.into(),
                    px_advance: item.px_advance as _,
                }
            }
        }

        pub struct FontAtlas {
            pub gpu: Arc<Gpu>,
            pub texture: Arc<Texture>,
            pub texture_data: Vec<u8>,
            pub texture_size: Vec2,
            pub map_char_info: HashMap<char, GlyphInfo>,
            pub px_scale: usize,
            pub cursor: (usize, usize),
            // pub dirty: bool,  // or even tracking the dirty areas directly for partial updates
            // most of the API for online rasterization is missing
            // COMPLEX USE CASES: multiple languages from multiple ab_fonts not necessarily all with the same scale (latin, chinese, arabic (rtl!), emojis)
            // potentially overflow mega-texture 8192x8192 (hardware limit) -> variable size allocator with LRU cache on a char render basis (but upstream data needed...))
        }

        impl FontAtlas {

            pub const BPP: u8 = 1;  // Bytes Per Pixel (equals number of channels for canonical u8 channels).
            pub const HOR_PADDING: u8 = 1;  // Horizontal padding between glyphs bounding box in the same row.
            pub const VER_PADDING: u8 = 1;  // Vertical padding between glyphs bounding box in consecutive rows. 
            pub const MAX_TEXTURE_DIM: u16 = 8192;  // Hardware texture dimension limit (relevant only for gpu textures).
    
            #[inline]
            pub fn stride(&self) -> usize {
                self.texture_size.x as usize * Self::BPP as usize
            }
    
            pub fn new<F: abFont>(gpu: Arc<Gpu>, ab_font: &F, px_scale: usize, chars: &[char]) -> Self {
                let scale = px_scale as f32;
                let ab_font_s = ab_font.as_scaled(scale);
                let ver_ascent = ab_font_s.ascent();

                let est_h = (scale + scale - ver_ascent).ceil();  // single glyph max height
                let tex_w = Self::MAX_TEXTURE_DIM as usize;
                let est_rows = (chars.len() as f32 / (tex_w as f32 / est_h)).ceil() as usize;  // using glyph max height (higher than width for glyphs) as proxy for max width (since we want upper bound anyways)
                let tex_h = (est_rows + 1) * (est_h as usize + Self::VER_PADDING as usize);  // height upper bound

                let num_channels = Self::BPP as usize;
                let stride = tex_w * num_channels;

                let mut texture_data = vec![0_u8; tex_h * stride];
                let mut map_char_info = HashMap::with_capacity(chars.len());
                let mut max_h_by_row = Vec::with_capacity(est_rows + 1);
                
                let mut cur_x = 1;
                let mut cur_y = 1;
                let mut cur_i = 0;
                max_h_by_row.push(0);
                
                // Pixel at (0, 0) reserved for white
                for channel_index in 0 .. num_channels { 
                    texture_data[channel_index] = 255;
                }
                
                for ch in chars.iter().cloned() {
                    let glyph_id = ab_font.glyph_id(ch);
                    let glyph = glyph_id.with_scale_and_position(scale, ab_glyph::point(0.0, 0.0));
                    let px_region = if let Some(out) = ab_font.outline_glyph(glyph) {
                        let bounds = out.px_bounds();
                        let cur_w = bounds.width().ceil() as usize;
                        let cur_h = bounds.height().ceil() as usize;
    
                        if (cur_x + cur_w) > tex_w {
                            cur_x = 0;
                            cur_y += max_h_by_row[cur_i] + Self::VER_PADDING as usize;
                            cur_i += 1;
                            max_h_by_row.push(0);
                        }
    
                        let ver_align = (ver_ascent + bounds.min.y).round().max(0.0) as usize;
                        let off_index = (cur_y + ver_align) as usize * stride + cur_x as usize;

                        out.draw(|x, y, c| {
                            let index = off_index + y as usize * stride + x as usize;
                            let value = (c * 255.0) as u8;
                            
                            for channel_index in 0 .. num_channels { 
                                texture_data[index + channel_index] = value;
                            }
                        });
    
                        let region = Rect::from_ltrb(
                            cur_x as _,
                            cur_y as _,
                            (cur_x + cur_w + 1) as _,
                            (cur_y + ver_align + cur_h + 1) as _,
                        );

                        let cur_h_from_last_row_bottom = ver_align + cur_h + 1;
                        if cur_h_from_last_row_bottom > max_h_by_row[cur_i] {
                            max_h_by_row[cur_i] = cur_h_from_last_row_bottom;
                        }
    
                        cur_x += cur_w + 1 + Self::HOR_PADDING as usize;
    
                        region
                    } else {
                        Rect::ZERO
                    };
                    let px_advance = ab_font_s.h_advance(glyph_id).ceil() as _;
    
                    map_char_info.insert(ch, GlyphInfo {
                        ab_id: glyph_id,
                        px_region,
                        px_advance,
                    });
                }

                let texture = match num_channels {
                    1 => Texture::from_a8(gpu.clone(), &texture_data, tex_w as _, tex_h as _),
                    4 => Texture::from_rgba8(gpu.clone(), &texture_data, tex_w as _, tex_h as _),
                    _ => unimplemented!(),
                };

                Self {
                    gpu,
                    texture: Arc::new(texture),
                    texture_data,
                    texture_size: (tex_w as _, tex_h as _).into(),
                    map_char_info,
                    px_scale,
                    cursor: (cur_x as _, cur_y as _).into()
                }
            }
            
        }
        
        pub struct Font {
            pub ab_font: ab_glyph::FontArc,
            pub monospace: Option<f32>,
            pub atlas: FontAtlas,
            pub char_fallback: char,
        }

        impl Font {

            pub const DEFAULT_CHAR_FALLBACK: char = '#';

            pub fn load(gpu: Arc<Gpu>, path: &str, monospace: bool) -> Self {
                let ascii = (32 ..= 126_u8).map(|c| c as char).collect::<Vec<_>>();
                let bytes = std::fs::read(path).unwrap();
                let ab_font = ab_glyph::FontArc::try_from_vec(bytes).unwrap();
    
                let mut atlas = FontAtlas::new(gpu, &ab_font, 32, &ascii);
                let space = atlas.map_char_info.get(&' ').unwrap().clone();

                assert!(space.px_region == Rect::ZERO, "space char incorrectly rasterized!");

                atlas.map_char_info.insert('\t', GlyphInfo {
                    ab_id: ab_font.glyph_id('\t'),
                    px_region: Rect::ZERO,
                    px_advance: space.px_advance * 4.0,
                });

                let char_fallback = Self::DEFAULT_CHAR_FALLBACK;
                if atlas.map_char_info.get(&char_fallback).is_none() {
                    todo!("atlas.rasterize(char_fallback)");
                }

                let monospace_adv = if monospace {
                    Some(space.px_advance)
                } else {
                    None
                };

                Self {
                    ab_font,
                    atlas,
                    monospace: monospace_adv,
                    char_fallback,
                }
            }
        
        }

        impl FontDrawInfo for Font {

            #[inline]
            fn px_monospace_advance(&self) -> Option<f32> {
                self.monospace
            }

            #[inline]
            fn px_kern(&self, pre_char: char, cur_char: char) -> f32 {
                if let Some(pre) = self.atlas.map_char_info.get(&pre_char) {
                    if let Some(cur) = self.atlas.map_char_info.get(&cur_char) {
                        return self.ab_font
                            .as_scaled(self.atlas.px_scale as f32)
                            .kern(pre.ab_id, cur.ab_id)
                    }
                }

                0.0
            }

            #[inline]
            fn texture_size(&self) -> Vec2 {
                self.atlas.texture_size
            }

            #[inline]
            fn get_char_draw_info_or_fallback(&mut self, ch: char) -> CharDrawInfo {
                match self.atlas.map_char_info.get(&ch) {
                    Some(char_info) => char_info.into(),
                    None => self.atlas.map_char_info.get(&self.char_fallback).unwrap().into(),
                }
            }
        }

    }

}


pub mod backend {
    use std::sync::Arc;
    use crate::core::draw::DrawOp;

    #[derive(Debug, Clone, Copy)]
    pub enum WindowMode {
        Windowed(f32, f32),
        Fullscreen,
    }
    
    impl Default for WindowMode {
        fn default() -> Self {
            Self::Windowed(800.0, 600.0)
        }
    }
        
    #[derive(Debug, Clone)]
    pub struct AppConfig {
        pub window_title: String,
        pub window_mode: WindowMode,
        pub window_resizable: bool,
        pub window_min_size: Option<(f32, f32)>,
        pub window_max_size: Option<(f32, f32)>,
        pub window_icon: Option<winit::window::Icon>,
        pub window_decorations: bool,
        pub cursor_icon: Option<winit::window::CursorIcon>,
        pub prefer_high_power: bool,
    }
    
    impl Default for AppConfig {
        fn default() -> Self {
            Self {
                window_title: String::from("Demo"),
                window_mode: Default::default(),
                window_resizable: true,
                window_min_size: Some((400.0, 300.0)),
                window_max_size: None,
                window_icon: None,
                window_decorations: true,
                cursor_icon: Some(Default::default()),
                prefer_high_power: false,
            }
        }
    }

    #[derive(Debug)]
    pub struct WinitWindowWrapper {
        pub raw: winit::window::Window,
        pub event_loop: Option<winit::event_loop::EventLoop<()>>,
    }
    
    #[derive(Debug)]
    pub struct WgpuSurfaceWrapper {
        pub raw: wgpu::Surface,
        pub config: wgpu::SurfaceConfiguration,
        pub clear_color: wgpu::Color,
    }
    
    #[derive(Debug)]
    pub struct Gpu {
        pub device: wgpu::Device,
        pub queue: wgpu::Queue,
    }
    
    #[derive(Debug, Default)]
    pub struct Renderer {
        pub draw_ops: Vec<DrawOp>,
    }

    impl Renderer {

        pub fn new() -> Self {
            Self {
                draw_ops: Vec::new(),
            }
        }

        pub fn push_draw_op(&mut self, draw_op: DrawOp) {
            self.draw_ops.push(draw_op);
        }
            
        pub fn submit_draw_ops<'slf: 'cmd_buf, 'rpass, 'cmd_buf>(&'slf mut self, render_pass: &'rpass mut wgpu::RenderPass<'cmd_buf>) {
            // TODO: optimize order (possibly already at push time) to minimize number of state changes.
            // i.e. same pipeline first, then same texture bind group, etc.

            for DrawOp {
                shader,
                uniform_buffers: _,
                texture_views_samplers: _,
                uniform_bind_group,
                texture_bind_group,
                vertex_buffer,
                index_buffer,
                draw_range
            } in self.draw_ops.iter() {

                // if pipeline != set_pipeline {
                    render_pass.set_pipeline(&shader.pipeline);
                //    set_pipeline = pipeline;
                // }

                let mut group_index = 0;

                // if uniform_bind_group != set_uniform_bind_group {
                if shader.uniform_bind_group_layout.is_some() {
                    render_pass.set_bind_group(group_index, uniform_bind_group.as_ref().unwrap(), &[]);
                    group_index += 1
                //    set_uniform_bind_group = uniform_bind_group;
                // }
                }

                // if texture_bind_group != set_texture_bind_group {
                if shader.texture_bind_group_layout.is_some() {
                    render_pass.set_bind_group(group_index, texture_bind_group.as_ref().unwrap(), &[]);
                //    set_texture_bind_group = texture_bind_group;
                // }
                }

                render_pass.set_vertex_buffer(0, vertex_buffer.as_ref().unwrap().slice(..));

                if let Some(index_buf) = &index_buffer {
                    render_pass.set_index_buffer(index_buf.slice(..), wgpu::IndexFormat::Uint16);
                    render_pass.draw_indexed(draw_range.start as _ .. draw_range.end as _, 0, 0 .. 1);
                } else {
                    render_pass.draw(draw_range.start as _ .. draw_range.end as _, 0 .. 1);
                }
            }
        }

        pub fn clear_draw_ops(&mut self) {
            self.draw_ops.clear();
        }

    }

    #[derive(Debug)]
    pub struct Backend {
        pub window: WinitWindowWrapper,
        pub surface: WgpuSurfaceWrapper,
        pub renderer: Renderer,
        pub gpu: Arc<Gpu>,
    }
    
    impl Backend {
    
        pub fn new(mut config: AppConfig) -> Self {
            env_logger::init();
        
            let instance = wgpu::Instance::new(wgpu::Backends::PRIMARY);
            let event_loop = winit::event_loop::EventLoop::new();
            
            #[allow(unused_mut)]
            let mut window_builder = winit::window::WindowBuilder::new()
                .with_title(&config.window_title)
                .with_window_icon(config.window_icon.take())
                .with_decorations(config.window_decorations)
                .with_resizable(config.window_resizable)
                .with_visible(true);
            
            match config.window_mode {
                WindowMode::Windowed(w, h) => {
                    window_builder = window_builder.with_inner_size(winit::dpi::LogicalSize::new(w, h));
                },
                WindowMode::Fullscreen => {
                    window_builder = window_builder.with_fullscreen(Some(winit::window::Fullscreen::Borderless(None)));
                },
            }
            
            if let Some((w, h)) = config.window_min_size {
                window_builder = window_builder.with_min_inner_size(winit::dpi::LogicalSize::new(w, h));
            }
            if let Some((w, h)) = config.window_max_size {
                window_builder = window_builder.with_max_inner_size(winit::dpi::LogicalSize::new(w, h));
            }
    
            let power_preference = if config.prefer_high_power { wgpu::PowerPreference::HighPerformance } else { wgpu::PowerPreference::LowPower };
    
            #[cfg(target_arch = "wasm32")]
            {
                use wasm_bindgen::JsCast;
                use winit::platform::web::WindowBuilderExtWebSys;
    
                let window = web_sys::window().unwrap();
                let document = window.document().unwrap();
                let canvas = document.query_selector("#canvas")
                    .expect("could not query for canvas!")
                    .expect("no #canvas found in document!")
                    .dyn_into::<web_sys::HtmlCanvasElement>()
                    .ok();
                
                let window_builder = window_builder.with_canvas(canvas);
            }
    
            #[allow(unused_mut)]
            let mut window = window_builder
                .build(&event_loop)
                .expect("could not build winit window!");
            
            if let Some(cursor_icon) = config.cursor_icon {
                window.set_cursor_visible(true);
                window.set_cursor_icon(cursor_icon);
            } else {
                window.set_cursor_visible(false);
            }
    
            #[cfg(target_arch = "wasm32")]
            {
                use winit::platform::web::WindowExtWebSys;
    
                let canvas = window.canvas();
    
                let window = web_sys::window().unwrap();
                let document = window.document().unwrap();
                let body = document.body().unwrap();
    
                body.append_child(&canvas).expect("could not append canvas to HTML body.");
            }
    
            let phys_size = window.inner_size();
    
            let surface = unsafe { instance.create_surface(&window) };
    
            let adapter = futures::executor::block_on(
                instance.request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference,
                    compatible_surface: Some(&surface),
                    ..Default::default()
                })
            ).unwrap();
            let (device, queue) = futures::executor::block_on(
                adapter.request_device(&wgpu::DeviceDescriptor {
                    label: None,
                    features: wgpu::Features::empty(),
                    limits: wgpu::Limits {
                        max_bind_groups: 2,
                        ..Default::default()
                    }
                }, None)
            ).unwrap();
    
            let surface_config = wgpu::SurfaceConfiguration {
                usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
                format: wgpu::TextureFormat::Bgra8UnormSrgb,
                width: phys_size.width as _,
                height: phys_size.height as _,
                present_mode: wgpu::PresentMode::Fifo,
            };
    
            surface.configure(&device, &surface_config);
    
            // Note: instance and adapter do not need to be kept alive.
            let gpu = Arc::new(Gpu {
                device,
                queue,
            });
    
            Self {
                window: WinitWindowWrapper {
                    raw: window,
                    event_loop: Some(event_loop),
                },
                surface: WgpuSurfaceWrapper {
                    raw: surface,
                    config: surface_config,
                    clear_color: wgpu::Color { r: 0.0, g: 0.0, b: 0.0, a: 0.0 },
                },
                renderer: Renderer::new(),
                gpu,
            }
        }
  
        #[inline]
        pub fn resize_surface(&mut self, size: winit::dpi::PhysicalSize<u32>) {
            self.surface.config.width = size.width;
            self.surface.config.height = size.height;
            self.surface.raw.configure(&self.gpu.device, &self.surface.config); 
        }

        #[inline]
        pub fn request_frame(&self) -> Option<wgpu::SurfaceTexture> {
            match self.surface.raw.get_current_texture() {
                Ok(frame) => Some(frame),
                Err(wgpu::SurfaceError::Lost) => {
                    self.surface.raw.configure(&self.gpu.device, &self.surface.config); 
                    match self.surface.raw.get_current_texture() {
                        Ok(frame) => Some(frame),
                        Err(e) => {
                            eprintln!("{:?}", e);
                            None
                        },
                    }
                },
                Err(e) => {
                    eprintln!("{:?}", e);
                    None
                },
            }
        }

    }

    impl Default for Backend {
        fn default() -> Self {
            Self::new(Default::default())
        }
    }

}


pub mod input {
    use winit::event::{WindowEvent, ElementState, MouseButton, MouseScrollDelta, TouchPhase};

    pub type Key = winit::event::VirtualKeyCode;
    pub type Mods = winit::event::ModifiersState;
    
    pub struct Input {
        pub mask: usize,
        pub mods: [Mods; 2],
        pub keys: [[i8; 256]; 2],
        pub mouse_buttons: [[i8; 256]; 2],
        pub mouse_position: [(f32, f32); 2],
        pub mouse_scroll_delta: f32,
        pub last_time_mouse_move: [std::time::Instant; 2],
        pub touch_state: [i8; 2], 
        pub touch_position: [(f32, f32); 2],
        pub last_time_touch_move: [std::time::Instant; 2],
        pub touch_moves_cursor_enabled: bool,
        pub chars_received: Option<Vec<char>>,
        // TODO: File drop
        // EXT: Higher level, drag? swipe? joystick? Allow customization?
    }
    
    impl Input {
    
        pub const MOUSE_SCROLL_STEP: f32 = 38.0;
        pub const MOUSE_BUTTON_LEFT: u8 = 0;
        pub const MOUSE_BUTTON_RIGHT: u8 = 1;
        pub const MOUSE_BUTTON_WHEEL: u8 = 2;
    
        pub fn new(touch_moves_cursor_enabled: bool) -> Self {
            use winit::event::ModifiersState;
    
            let now = std::time::Instant::now();
    
            Self {
                mask: 0,
                mods: [ModifiersState::empty(), ModifiersState::empty()],
                keys: [[0; 256], [0; 256]],
                mouse_buttons: [[0; 256], [0; 256]],
                mouse_position: [(0.0, 0.0), (0.0, 0.0)],
                mouse_scroll_delta: 0.0,
                last_time_mouse_move: [now, now],
                touch_state: [0, 0],
                touch_position: [(0.0, 0.0), (0.0, 0.0)],
                last_time_touch_move: [now, now],
                touch_moves_cursor_enabled,
                chars_received: Some(Vec::new()),
            }
        }
    
        // modifier keys
        #[inline]
        pub fn ctrl_down(&self) -> bool {
            self.mods[self.mask].ctrl()
        }
    
        #[inline]
        pub fn alt_down(&self) -> bool {
            self.mods[self.mask].alt()
        }
    
        #[inline]
        pub fn shift_down(&self) -> bool {
            self.mods[self.mask].shift()
        }
    
        #[inline]
        pub fn logo_down(&self) -> bool {
            self.mods[self.mask].logo()
        }
    
        #[inline]
        pub fn ctrl_up(&self) -> bool {
            !self.mods[self.mask].ctrl()
        }
    
        #[inline]
        pub fn alt_up(&self) -> bool {
            !self.mods[self.mask].alt()
        }
    
        #[inline]
        pub fn shift_up(&self) -> bool {
            !self.mods[self.mask].shift()
        }
    
        #[inline]
        pub fn logo_up(&self) -> bool {
            !self.mods[self.mask].logo()
        }
    
        #[inline]
        pub fn ctrl_just_down(&self) -> bool {
            self.mods[self.mask].ctrl() && !self.mods[1 - self.mask].ctrl()
        }
    
        #[inline]
        pub fn alt_just_down(&self) -> bool {
            self.mods[self.mask].alt() && !self.mods[1 - self.mask].alt()
        }
    
        #[inline]
        pub fn shift_just_down(&self) -> bool {
            self.mods[self.mask].shift() && !self.mods[1 - self.mask].shift()
        }
    
        #[inline]
        pub fn logo_just_down(&self) -> bool {
            self.mods[self.mask].logo() && !self.mods[1 - self.mask].logo()
        }
    
        #[inline]
        pub fn ctrl_just_up(&self) -> bool {
            !self.mods[self.mask].ctrl() && self.mods[1 - self.mask].ctrl()
        }
    
        #[inline]
        pub fn alt_just_up(&self) -> bool {
            !self.mods[self.mask].alt() && self.mods[1 - self.mask].alt()
        }
    
        #[inline]
        pub fn shift_just_up(&self) -> bool {
            !self.mods[self.mask].shift() && self.mods[1 - self.mask].shift()
        }
    
        #[inline]
        pub fn logo_just_up(&self) -> bool {
            !self.mods[self.mask].logo() && self.mods[1 - self.mask].logo()
        }
    
        // normal keys
        #[inline]
        pub fn key_down(&self, key: Key) -> bool {
            let key_index = key as usize;
            self.keys[self.mask][key_index] > 0
        }
    
        #[inline]
        pub fn key_up(&self, key: Key) -> bool {
            let key_index = key as usize;
            self.keys[self.mask][key_index] < 1
        }
    
        #[inline]
        pub fn key_just_down(&self, key: Key) -> bool {
            let key_index = key as usize;
            (self.keys[self.mask][key_index] - self.keys[1 - self.mask][key_index]) > 0
        }
    
        #[inline]
        pub fn key_just_up(&self, key: Key) -> bool {
            let key_index = key as usize;
            (self.keys[self.mask][key_index] - self.keys[1 - self.mask][key_index]) < 0
        }
    
        // mouse
        #[inline]
        pub fn mouse_down(&self, button: u8) -> bool {
            self.mouse_buttons[self.mask][button as usize] > 0
        }
    
        #[inline]
        pub fn mouse_up(&self, button: u8) -> bool {
            self.mouse_buttons[self.mask][button as usize] < 1
        }
    
        #[inline]
        pub fn mouse_just_down(&self, button: u8) -> bool {
            (self.mouse_buttons[self.mask][button as usize] - self.mouse_buttons[1 - self.mask][button as usize]) > 0
        }
    
        #[inline]
        pub fn mouse_just_up(&self, button: u8) -> bool {
            (self.mouse_buttons[self.mask][button as usize] - self.mouse_buttons[1 - self.mask][button as usize]) < 0
        }
    
        #[inline]
        pub fn mouse_pos(&self) -> (f32, f32) {
            self.mouse_position[self.mask]
        }
    
        #[inline]
        pub fn prev_mouse_pos(&self) -> (f32, f32) {
            self.mouse_position[1 - self.mask]
        }
    
        #[inline]
        pub fn mouse_delta(&self) -> (f32, f32) {
            let (prev_x, prev_y) = self.mouse_position[1 - self.mask];
            let (curr_x, curr_y) = self.mouse_position[self.mask];
            (curr_x - prev_x, curr_y - prev_y)
        }
    
        #[inline]
        pub fn mouse_scroll(&self) -> f32 {
            self.mouse_scroll_delta
        }
    
        #[inline]
        pub fn get_mouse_just_down(&self, button: u8) -> Option<(f32, f32)> {
            if self.mouse_just_down(button) {
                Some(self.mouse_pos())
            } else {
                None
            }
        }
    
        #[inline]
        pub fn get_mouse_just_up(&self, button: u8) -> Option<(f32, f32)> {
            if self.mouse_just_up(button) {
                Some(self.mouse_pos())
            } else {
                None
            }
        }
    
        // touch
        #[inline]
        pub fn touch_down(&self) -> bool {
            self.touch_state[self.mask] > 0
        }
    
        #[inline]
        pub fn touch_up(&self) -> bool {
            self.touch_state[self.mask] < 1
        }
    
        #[inline]
        pub fn touch_just_down(&self) -> bool {
            self.touch_state[self.mask] > 0 && self.touch_state[1 - self.mask] < 1
        }
    
        #[inline]
        pub fn touch_just_up(&self) -> bool {
            self.touch_state[self.mask] < 1 && self.touch_state[1 - self.mask] > 0
        }
    
        #[inline]
        pub fn touch_pos(&self) -> (f32, f32) {
            self.touch_position[self.mask]
        }
    
        #[inline]
        pub fn prev_touch_pos(&self) -> (f32, f32) {
            self.touch_position[1 - self.mask]
        }
    
        #[inline]
        pub fn touch_delta(&self) -> (f32, f32) {
            let (prev_x, prev_y) = self.touch_position[1 - self.mask];
            let (curr_x, curr_y) = self.touch_position[self.mask];
            (curr_x - prev_x, curr_y - prev_y)
        }
    
        #[inline]
        pub fn get_touch_just_down(&self) -> Option<(f32, f32)> {
            if self.touch_just_down() {
                Some(self.touch_pos())
            } else {
                None
            }
        }
    
        #[inline]
        pub fn get_touch_just_up(&self) -> Option<(f32, f32)> {
            if self.touch_just_up() {
                Some(self.touch_pos())
            } else {
                None
            }
        }
    
        // hybrid mouse and touch together
        #[inline]
        pub fn user_down(&self) -> bool {
            self.mouse_down(Self::MOUSE_BUTTON_LEFT) || self.touch_down()
        }
    
        #[inline]
        pub fn user_up(&self) -> bool {
            self.mouse_up(Self::MOUSE_BUTTON_LEFT) || self.touch_up()
        }
    
        #[inline]
        pub fn user_just_down(&self) -> bool {
            self.mouse_just_down(Self::MOUSE_BUTTON_LEFT) || self.touch_just_down()
        }
    
        #[inline]
        pub fn user_just_up(&self) -> bool {
            self.mouse_just_up(Self::MOUSE_BUTTON_LEFT) || self.touch_just_up()
        }
    
        #[inline]
        pub fn user_pos(&self) -> (f32, f32) {
            if self.last_time_touch_move[self.mask] > self.last_time_mouse_move[self.mask] {
                self.touch_pos()
            } else {
                self.mouse_pos()
            }
        }
    
        #[inline]
        pub fn prev_user_pos(&self) -> (f32, f32) {
            if self.last_time_touch_move[1 - self.mask] > self.last_time_mouse_move[1 - self.mask] {
                self.touch_pos()
            } else {
                self.mouse_pos()
            }
        }
    
        #[inline]
        pub fn user_delta(&self) -> (f32, f32) {
            let (prev_x, prev_y) = self.prev_user_pos();
            let (curr_x, curr_y) = self.user_pos();
            (curr_x - prev_x, curr_y - prev_y)
        }
    
        #[inline]
        pub fn get_user_just_down(&self) -> Option<(f32, f32)> {
            let mouse_at = self.get_mouse_just_down(Self::MOUSE_BUTTON_LEFT);
    
            if mouse_at.is_none() {
                self.get_touch_just_down()
            } else {
                mouse_at
            }
        }
    
        #[inline]
        pub fn get_user_just_up(&self) -> Option<(f32, f32)> {
            let mouse_at = self.get_mouse_just_up(Self::MOUSE_BUTTON_LEFT);
    
            if mouse_at.is_none() {
                self.get_touch_just_up()
            } else {
                mouse_at
            }
        }
    
        // key input
        #[inline]
        pub fn take_chars_received(&mut self) -> Option<Vec<char>> {
            self.chars_received.take()
        }
    
        // swap buffers and reset mouse_scroll_delta and chars_received
        #[inline]
        pub fn update_frame(&mut self) {
            // swap buffers (writes now affect at [mask])
            self.mask = 1 - self.mask;
    
            // reset mouse scroll (which is on a per frame basis, accumulation should happen upstream as needed)
            self.mouse_scroll_delta = 0.0;
    
            // reset chars-received buffer
            self.chars_received = Some(Vec::new());
        }
    
        // handle event and update state
        #[inline]
        pub fn handle_window_event(&mut self, window: &mut winit::window::Window, event: winit::event::WindowEvent) {
            match event {
                WindowEvent::ModifiersChanged(mods) => {
                    self.mods[self.mask] = mods;
                },
                WindowEvent::KeyboardInput { ref input, .. } => {
                    if let Some(key_code) = input.virtual_keycode {
                        let key_index = key_code as usize;
                        self.keys[self.mask][key_index] = match input.state {
                            ElementState::Pressed => 1,
                            ElementState::Released => 0,
                        };
                    }
                },
                WindowEvent::CursorMoved { position, .. } => {
                    self.mouse_position[self.mask] = position.to_logical::<f32>(window.scale_factor()).into();
                    self.last_time_mouse_move[self.mask] = std::time::Instant::now();
                },
                WindowEvent::Touch(touch) => {
                    let abs_force: i8 = {
                        if let Some(force) = touch.force {
                            (force.normalized() * i8::MAX as f64).round() as i8
                        } else {
                            1
                        }
                    };
    
                    match touch.phase {
                        TouchPhase::Started => {
                            self.touch_state[self.mask] = abs_force;
                        },
                        TouchPhase::Ended | TouchPhase::Cancelled => {
                            self.touch_state[self.mask] = -abs_force;
                        },
                        TouchPhase::Moved => {
                            self.touch_state[self.mask] = (self.touch_state[self.mask]).signum() * abs_force;
                        },
                    };
    
                    self.touch_position[self.mask] = touch.location.to_logical::<f32>(window.scale_factor()).into();
                    self.last_time_touch_move[self.mask] = std::time::Instant::now();
                    if self.touch_moves_cursor_enabled {
                        if let Err(e) = window.set_cursor_position(touch.location) {
                            eprintln!("could not set cursor position to touch location due to {:?}", e);
                        }
                    }
                },
                WindowEvent::CursorLeft { .. } => {
                    self.mouse_position[self.mask] = (-1.0, -1.0);
                },
                WindowEvent::MouseInput { state, button, .. } => {
                    let s = match state {
                        ElementState::Pressed => 1,
                        ElementState::Released => 0,
                    };
    
                    match button {
                        MouseButton::Left => self.mouse_buttons[self.mask][0] = s,
                        MouseButton::Right => self.mouse_buttons[self.mask][1] = s,
                        MouseButton::Middle => self.mouse_buttons[self.mask][2] = s,
                        MouseButton::Other(byte) => self.mouse_buttons[self.mask][byte as u8 as usize] = s,
                    };
    
                    self.last_time_mouse_move[self.mask] = std::time::Instant::now();
                },
                WindowEvent::MouseWheel { delta, .. } => {
                    match delta {
                        MouseScrollDelta::LineDelta(_, y) => {
                            self.mouse_scroll_delta += y;
                        },
                        MouseScrollDelta::PixelDelta(delta) => {
                            self.mouse_scroll_delta += (delta.y as f32 / Self::MOUSE_SCROLL_STEP) as f32;
                        },
                    };
                },
                WindowEvent::ReceivedCharacter(c) => {
                    if let Some(vec_chars) = &mut self.chars_received {
                        vec_chars.push(c);
                    } else {
                        panic!("update_frame() not called or take_chars_received() called during event handling (wait for all chars to be buffered!)");
                    }
                },
                WindowEvent::DroppedFile(_path_buf) => {
                    todo!()
                },
                WindowEvent::HoveredFile(_path_buf) => {
                    todo!()
                },
                WindowEvent::HoveredFileCancelled => {
                    todo!()
                },
                _ => {
    
                },
            };
        }
    
    }
    
    impl Default for Input {
        fn default() -> Self {
            Self::new(true)
        }
    }
    
}


pub mod clock {

    #[derive(Debug)]
    pub struct Clock {
        pub last_frame_instant: std::time::Instant,
        pub last_frame_duration: std::time::Duration,
    }

    impl Clock {

        #[inline]
        pub fn update_frame(&mut self) {
            let cur_frame_instant = std::time::Instant::now();
            self.last_frame_duration = cur_frame_instant - self.last_frame_instant;
            self.last_frame_instant = cur_frame_instant;
        }

        #[inline]
        pub fn delta_time(&mut self) -> std::time::Duration {
            self.last_frame_duration
        }

    }

    impl Default for Clock {
        fn default() -> Self {
            Self {
                last_frame_instant: std::time::Instant::now(),
                last_frame_duration: std::time::Duration::new(0, 0),
            }
        }
    }

}


pub mod app {
    use ::core::any::TypeId;
    use crate::core::any_map::AnyMap;
    use crate::backend::*;
    use crate::input::Input;
    use crate::clock::Clock;
    use std::sync::Arc;

    pub struct App {
        pub any_map: AnyMap,
        pub backend: Backend,
        pub clock: Clock,
        pub input: Input,
        pub update: Option<Box<dyn FnMut(&mut Self)>>,
        pub render: Option<Box<dyn FnMut(&mut Self)>>,
        pub exited: bool,
    }

    impl App {

        pub fn new(config: AppConfig) -> Self {
            Self {
                any_map: Default::default(),
                backend: Backend::new(config),
                clock: Default::default(),
                input: Default::default(),
                update: None,
                render: None,
                exited: false,
            }
        }

        #[inline]
        pub fn insert<T: 'static>(&mut self, resource: T) -> Option<T> {
            if let Some(res) = self.any_map.insert(TypeId::of::<T>(), Box::new(resource)) {
                Some(*res.downcast::<T>().unwrap())
            } else {
                None
            }
        }

        #[inline]
        pub fn remove<T: 'static>(&mut self) -> Option<T> {
            if let Some(res) = self.any_map.remove(&TypeId::of::<T>()) {
                Some(*res.downcast::<T>().unwrap())
            } else {
                None
            }
        }

        #[inline]
        pub fn get<T: 'static>(&self) -> Option<&T> {
            if let Some(res) = self.any_map.get(&TypeId::of::<T>()) {
                res.downcast_ref::<T>()
            } else {
                None
            }
        }

        #[inline]
        pub fn get_mut<T: 'static>(&mut self) -> Option<&mut T> {
            if let Some(res) = self.any_map.get_mut(&TypeId::of::<T>()) {
                res.downcast_mut::<T>()   
            } else {
                None
            }
        }
        
        pub fn update<F: FnMut(&mut Self) + 'static>(&mut self, function: F) {
            self.update = Some(Box::new(function));        
        }

        pub fn render<F: FnMut(&mut Self) + 'static>(&mut self, function: F) {
            self.render = Some(Box::new(function));
        }

        #[inline]
        pub fn gpu(&self) -> Arc<Gpu> {
            self.backend.gpu.clone()
        }

        #[inline]
        pub fn exit(&mut self) {
            self.exited = true;
        }

        #[inline]
        pub fn get_viewport_size(&self) -> (f32, f32) {
            let size = self.backend.window.raw.inner_size().to_logical(self.backend.window.raw.scale_factor());
            (size.width, size.height)
        }

        #[inline]
        pub fn set_title(&self, title: &str) {
            self.backend.window.raw.set_title(title);
        }

        #[inline]
        pub fn handle_window_event(&mut self, event: winit::event::WindowEvent) {
            match event {
                winit::event::WindowEvent::CloseRequested | winit::event::WindowEvent::Destroyed => {
                    self.exit();
                },
                winit::event::WindowEvent::Resized(size) => {
                    self.backend.resize_surface(size);
                },
                winit::event::WindowEvent::ScaleFactorChanged { scale_factor: _, new_inner_size } => {
                    self.backend.resize_surface(*new_inner_size);
                },
                _ => {
                    self.input.handle_window_event(&mut self.backend.window.raw, event);
                }
            };
        }

        #[inline]
        pub fn update_frame(&mut self) {
            self.clock.update_frame();
            self.input.update_frame();
            self.backend.window.raw.request_redraw();
        }

        #[inline]
        pub fn render_frame(&mut self) {
            if let Some(frame) = self.backend.request_frame() {
                let frame_view = frame.texture.create_view(&wgpu::TextureViewDescriptor::default());
                let mut command_encoder = self.backend.gpu.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
                let mut render_pass = command_encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: None,
                    color_attachments: &[
                        wgpu::RenderPassColorAttachment {
                            view: &frame_view,
                            resolve_target: None,
                            ops: wgpu::Operations {
                                load: wgpu::LoadOp::Clear(self.backend.surface.clear_color),
                                store: true,
                            },
                        },
                    ],
                    depth_stencil_attachment: None,
                });
                self.backend.renderer.submit_draw_ops(&mut render_pass);
                drop(render_pass);
                self.backend.gpu.queue.submit(std::iter::once(command_encoder.finish()));
                frame.present();
                self.backend.renderer.clear_draw_ops();
            } else {
                eprintln!("Failed to acquire next swap chain frame!");
                self.exit()
            }
        }

    }

    impl Default for App {
        fn default() -> Self {
            Self::new(Default::default())
        }
    }

    pub fn run(mut app: App) {
        use winit::event::Event;
        use winit::event_loop::{ControlFlow, EventLoopWindowTarget};

        let event_loop = app.backend.window.event_loop.take().unwrap();

        let mut update = app.update.take().unwrap_or(Box::new(|_app| {}));
        let mut render = app.render.take().unwrap_or(Box::new(|_app| {}));

        let main = move |event: Event<()>, _: &EventLoopWindowTarget<()>, control_flow: &mut ControlFlow| {
            *control_flow = ControlFlow::Poll;
        
            match event {
                Event::WindowEvent { event, window_id: _, .. } => {
                    app.handle_window_event(event);
                },
                Event::MainEventsCleared => {
                    app.update_frame();
                    update(&mut app);
    
                    if app.exited {
                        *control_flow = ControlFlow::Exit;
                    }
                },
                Event::RedrawRequested(_) => {
                    render(&mut app);
                    app.render_frame();
                },
                _ => (),
            
            };
        };

        event_loop.run(main);
    }

}


pub mod prelude {
    pub use crate::core::{math::*, color::*, draw::*, draw::prefabs::*};
    pub use crate::auxil::{*, draw::*, font::*};
    pub use crate::backend::{WindowMode, AppConfig};
    pub use crate::input::*;
    pub use crate::clock::*;
    pub use crate::app::*;
}
