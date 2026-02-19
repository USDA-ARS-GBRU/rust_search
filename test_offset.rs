use std::os::raw::{c_char, c_double, c_int};

#[repr(C)]
pub struct ThalResultsC {
    pub msg: [c_char; 255],
    pub temp: c_double,
}

fn main() {
    println!("{}", std::mem::offset_of!(ThalResultsC, temp));
}
