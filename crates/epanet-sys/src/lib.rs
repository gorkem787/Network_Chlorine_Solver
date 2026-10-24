//! Raw bindings to the subset of the EPANET 2.2 thread-safe (`EN_`) API that
//! this project uses. Signatures follow `EPANET/include/epanet2_2.h`; constant
//! values follow `epanet2_enums.h`.

#![allow(non_camel_case_types)]

use std::ffi::{c_char, c_double, c_int, c_long, c_void};

pub type EN_Project = *mut c_void;

pub const EN_MAXID: usize = 31;
pub const EN_MAXMSG: usize = 255;

// Count types
pub const EN_NODECOUNT: c_int = 0;
pub const EN_LINKCOUNT: c_int = 2;

// Node properties
pub const EN_BASEDEMAND: c_int = 1;
pub const EN_INITQUAL: c_int = 4;
pub const EN_SOURCEQUAL: c_int = 5;
pub const EN_DEMAND: c_int = 9;
pub const EN_QUALITY: c_int = 12;
pub const EN_TANK_KBULK: c_int = 23;

// Link properties
pub const EN_DIAMETER: c_int = 0;
pub const EN_LENGTH: c_int = 1;
pub const EN_KBULK: c_int = 6;
pub const EN_KWALL: c_int = 7;
pub const EN_FLOW: c_int = 8;

// Time parameters
pub const EN_DURATION: c_int = 0;
pub const EN_HYDSTEP: c_int = 1;
pub const EN_QUALSTEP: c_int = 2;
pub const EN_REPORTSTEP: c_int = 5;
pub const EN_REPORTSTART: c_int = 6;

// Node types
pub const EN_JUNCTION: c_int = 0;
pub const EN_RESERVOIR: c_int = 1;
pub const EN_TANK: c_int = 2;

// Link types
pub const EN_CVPIPE: c_int = 0;
pub const EN_PIPE: c_int = 1;
pub const EN_PUMP: c_int = 2;

// Quality types
pub const EN_NONE: c_int = 0;
pub const EN_CHEM: c_int = 1;

// Flow units
pub const EN_CFS: c_int = 0;
pub const EN_GPM: c_int = 1;
pub const EN_MGD: c_int = 2;
pub const EN_IMGD: c_int = 3;
pub const EN_AFD: c_int = 4;
pub const EN_LPS: c_int = 5;
pub const EN_LPM: c_int = 6;
pub const EN_MLD: c_int = 7;
pub const EN_CMH: c_int = 8;
pub const EN_CMD: c_int = 9;

pub const EN_NOSAVE: c_int = 0;

extern "C" {
    pub fn EN_createproject(ph: *mut EN_Project) -> c_int;
    pub fn EN_deleteproject(ph: EN_Project) -> c_int;
    pub fn EN_open(ph: EN_Project, inp: *const c_char, rpt: *const c_char, out: *const c_char) -> c_int;
    pub fn EN_close(ph: EN_Project) -> c_int;
    pub fn EN_geterror(errcode: c_int, out_errmsg: *mut c_char, max_len: c_int) -> c_int;

    pub fn EN_getcount(ph: EN_Project, object: c_int, count: *mut c_int) -> c_int;
    pub fn EN_getflowunits(ph: EN_Project, units: *mut c_int) -> c_int;
    pub fn EN_gettimeparam(ph: EN_Project, param: c_int, value: *mut c_long) -> c_int;
    pub fn EN_settimeparam(ph: EN_Project, param: c_int, value: c_long) -> c_int;
    pub fn EN_getqualtype(ph: EN_Project, qual_type: *mut c_int, trace_node: *mut c_int) -> c_int;
    pub fn EN_setqualtype(
        ph: EN_Project,
        qual_type: c_int,
        chem_name: *const c_char,
        chem_units: *const c_char,
        trace_node: *const c_char,
    ) -> c_int;

    pub fn EN_getnodeid(ph: EN_Project, index: c_int, out_id: *mut c_char) -> c_int;
    pub fn EN_getnodetype(ph: EN_Project, index: c_int, node_type: *mut c_int) -> c_int;
    pub fn EN_getnodevalue(ph: EN_Project, index: c_int, property: c_int, value: *mut c_double) -> c_int;
    pub fn EN_setnodevalue(ph: EN_Project, index: c_int, property: c_int, value: c_double) -> c_int;
    pub fn EN_getcoord(ph: EN_Project, index: c_int, x: *mut c_double, y: *mut c_double) -> c_int;

    pub fn EN_getlinkid(ph: EN_Project, index: c_int, out_id: *mut c_char) -> c_int;
    pub fn EN_getlinktype(ph: EN_Project, index: c_int, link_type: *mut c_int) -> c_int;
    pub fn EN_getlinknodes(ph: EN_Project, index: c_int, node1: *mut c_int, node2: *mut c_int) -> c_int;
    pub fn EN_getlinkvalue(ph: EN_Project, index: c_int, property: c_int, value: *mut c_double) -> c_int;
    pub fn EN_setlinkvalue(ph: EN_Project, index: c_int, property: c_int, value: c_double) -> c_int;

    pub fn EN_openH(ph: EN_Project) -> c_int;
    pub fn EN_initH(ph: EN_Project, init_flag: c_int) -> c_int;
    pub fn EN_runH(ph: EN_Project, current_time: *mut c_long) -> c_int;
    pub fn EN_nextH(ph: EN_Project, t_step: *mut c_long) -> c_int;
    pub fn EN_closeH(ph: EN_Project) -> c_int;

    pub fn EN_openQ(ph: EN_Project) -> c_int;
    pub fn EN_initQ(ph: EN_Project, save_flag: c_int) -> c_int;
    pub fn EN_runQ(ph: EN_Project, current_time: *mut c_long) -> c_int;
    pub fn EN_nextQ(ph: EN_Project, t_step: *mut c_long) -> c_int;
    pub fn EN_closeQ(ph: EN_Project) -> c_int;
}
