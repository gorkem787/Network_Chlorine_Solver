//! Minimal safe wrapper around an EPANET project handle.

use std::ffi::{c_char, c_int, c_long, CStr, CString};
use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use epanet_sys as sys;

#[derive(Debug, Clone, PartialEq)]
pub struct EpanetError {
    /// EPANET error code, or 0 for errors raised by this crate.
    pub code: i32,
    pub message: String,
}

impl EpanetError {
    pub(crate) fn other(message: impl Into<String>) -> Self {
        EpanetError { code: 0, message: message.into() }
    }

    fn from_code(code: c_int, context: &str) -> Self {
        let mut buf = [0 as c_char; sys::EN_MAXMSG + 1];
        // SAFETY: buf holds EN_MAXMSG + 1 bytes and EPANET writes at most max_len.
        unsafe { sys::EN_geterror(code, buf.as_mut_ptr(), sys::EN_MAXMSG as c_int) };
        let msg = unsafe { CStr::from_ptr(buf.as_ptr()) }.to_string_lossy().into_owned();
        EpanetError { code, message: format!("{context}: {msg}") }
    }
}

impl fmt::Display for EpanetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.code > 0 {
            write!(f, "EPANET error {}: {}", self.code, self.message)
        } else {
            f.write_str(&self.message)
        }
    }
}

impl std::error::Error for EpanetError {}

/// EPANET codes below 100 are warnings (e.g. negative pressures) and the run
/// continues; 100 and above are errors.
fn check(code: c_int, context: &str) -> Result<(), EpanetError> {
    if code >= 100 {
        Err(EpanetError::from_code(code, context))
    } else {
        Ok(())
    }
}

pub struct Project {
    ph: sys::EN_Project,
    temp_files: Vec<PathBuf>,
}

static COUNTER: AtomicU64 = AtomicU64::new(0);

fn temp_path(ext: &str) -> PathBuf {
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!("chlorine_{}_{n}.{ext}", std::process::id()))
}

fn c_path(p: &Path) -> Result<CString, EpanetError> {
    CString::new(p.to_string_lossy().as_bytes()).map_err(|_| EpanetError::other("path contains a NUL byte"))
}

impl Project {
    pub fn open(inp: &Path) -> Result<Self, EpanetError> {
        // EPANET opens files with the C runtime's narrow fopen and limits paths
        // to 259 characters, so non-ASCII names (e.g. "hızlı_net.inp") and deep
        // folders fail. Work on a temp copy with a plain name instead.
        let inp_copy = temp_path("inp");
        std::fs::copy(inp, &inp_copy)
            .map_err(|e| EpanetError::other(format!("cannot read {}: {e}", inp.display())))?;
        let rpt = temp_path("rpt");

        let mut ph: sys::EN_Project = std::ptr::null_mut();
        // SAFETY: ph is a valid out-pointer.
        check(unsafe { sys::EN_createproject(&mut ph) }, "create project")?;
        let project = Project { ph, temp_files: vec![inp_copy.clone(), rpt.clone()] };

        let (inp_c, rpt_c, out_c) = (c_path(&inp_copy)?, c_path(&rpt)?, CString::default());
        // SAFETY: all strings are NUL-terminated and outlive the call.
        let code = unsafe { sys::EN_open(project.ph, inp_c.as_ptr(), rpt_c.as_ptr(), out_c.as_ptr()) };
        check(code, &format!("open {}", inp.display()))?;
        Ok(project)
    }

    pub fn count(&self, object: c_int) -> Result<usize, EpanetError> {
        let mut n = 0;
        check(unsafe { sys::EN_getcount(self.ph, object, &mut n) }, "get count")?;
        Ok(n as usize)
    }

    pub fn flow_units(&self) -> Result<c_int, EpanetError> {
        let mut u = 0;
        check(unsafe { sys::EN_getflowunits(self.ph, &mut u) }, "get flow units")?;
        Ok(u)
    }

    pub fn time_param(&self, param: c_int) -> Result<i64, EpanetError> {
        let mut v: c_long = 0;
        check(unsafe { sys::EN_gettimeparam(self.ph, param, &mut v) }, "get time parameter")?;
        Ok(v as i64)
    }

    pub fn set_time_param(&mut self, param: c_int, value: i64) -> Result<(), EpanetError> {
        let v = c_long::try_from(value).map_err(|_| EpanetError::other(format!("time value {value} s is too large")))?;
        check(unsafe { sys::EN_settimeparam(self.ph, param, v) }, "set time parameter")
    }

    pub fn qual_type(&self) -> Result<c_int, EpanetError> {
        let (mut q, mut trace) = (0, 0);
        check(unsafe { sys::EN_getqualtype(self.ph, &mut q, &mut trace) }, "get quality type")?;
        Ok(q)
    }

    pub fn set_chemical(&mut self, name: &str, units: &str) -> Result<(), EpanetError> {
        let (n, u, t) = (CString::new(name).unwrap(), CString::new(units).unwrap(), CString::default());
        check(
            unsafe { sys::EN_setqualtype(self.ph, sys::EN_CHEM, n.as_ptr(), u.as_ptr(), t.as_ptr()) },
            "set quality type",
        )
    }

    fn id(&self, index: usize, f: unsafe extern "C" fn(sys::EN_Project, c_int, *mut c_char) -> c_int) -> Result<String, EpanetError> {
        let mut buf = [0 as c_char; sys::EN_MAXID + 1];
        // SAFETY: EPANET writes at most EN_MAXID characters plus NUL.
        check(unsafe { f(self.ph, index as c_int + 1, buf.as_mut_ptr()) }, "get id")?;
        Ok(unsafe { CStr::from_ptr(buf.as_ptr()) }.to_string_lossy().into_owned())
    }

    // Node and link indices are 0-based here and 1-based in EPANET.

    pub fn node_id(&self, i: usize) -> Result<String, EpanetError> {
        self.id(i, sys::EN_getnodeid)
    }

    pub fn node_type(&self, i: usize) -> Result<c_int, EpanetError> {
        let mut t = 0;
        check(unsafe { sys::EN_getnodetype(self.ph, i as c_int + 1, &mut t) }, "get node type")?;
        Ok(t)
    }

    pub fn node_value(&self, i: usize, prop: c_int) -> Result<f64, EpanetError> {
        let mut v = 0.0;
        check(unsafe { sys::EN_getnodevalue(self.ph, i as c_int + 1, prop, &mut v) }, "get node value")?;
        Ok(v)
    }

    pub fn set_node_value(&mut self, i: usize, prop: c_int, v: f64) -> Result<(), EpanetError> {
        check(unsafe { sys::EN_setnodevalue(self.ph, i as c_int + 1, prop, v) }, "set node value")
    }

    /// `None` when the input file has no [COORDINATES] entry for the node.
    pub fn node_coords(&self, i: usize) -> Option<(f64, f64)> {
        let (mut x, mut y) = (0.0, 0.0);
        let code = unsafe { sys::EN_getcoord(self.ph, i as c_int + 1, &mut x, &mut y) };
        (code == 0).then_some((x, y))
    }

    pub fn link_id(&self, i: usize) -> Result<String, EpanetError> {
        self.id(i, sys::EN_getlinkid)
    }

    pub fn link_type(&self, i: usize) -> Result<c_int, EpanetError> {
        let mut t = 0;
        check(unsafe { sys::EN_getlinktype(self.ph, i as c_int + 1, &mut t) }, "get link type")?;
        Ok(t)
    }

    /// Start and end node, 0-based.
    pub fn link_nodes(&self, i: usize) -> Result<(usize, usize), EpanetError> {
        let (mut a, mut b) = (0, 0);
        check(unsafe { sys::EN_getlinknodes(self.ph, i as c_int + 1, &mut a, &mut b) }, "get link nodes")?;
        Ok((a as usize - 1, b as usize - 1))
    }

    pub fn link_value(&self, i: usize, prop: c_int) -> Result<f64, EpanetError> {
        let mut v = 0.0;
        check(unsafe { sys::EN_getlinkvalue(self.ph, i as c_int + 1, prop, &mut v) }, "get link value")?;
        Ok(v)
    }

    pub fn set_link_value(&mut self, i: usize, prop: c_int, v: f64) -> Result<(), EpanetError> {
        check(unsafe { sys::EN_setlinkvalue(self.ph, i as c_int + 1, prop, v) }, "set link value")
    }

    /// Runs hydraulics and water quality together, calling `record(t)` with
    /// EPANET's clock (seconds) after each hydraulic event has been solved.
    pub fn run(&mut self, mut record: impl FnMut(&Self, i64) -> Result<(), EpanetError>) -> Result<(), EpanetError> {
        unsafe {
            check(sys::EN_openH(self.ph), "open hydraulics")?;
            check(sys::EN_initH(self.ph, sys::EN_NOSAVE), "initialize hydraulics")?;
            check(sys::EN_openQ(self.ph), "open quality")?;
            check(sys::EN_initQ(self.ph, sys::EN_NOSAVE), "initialize quality")?;
        }
        let result = (|| loop {
            let (mut t, mut tq, mut step, mut qstep): (c_long, c_long, c_long, c_long) = (0, 0, 0, 0);
            unsafe {
                check(sys::EN_runH(self.ph, &mut t), "run hydraulics")?;
                check(sys::EN_runQ(self.ph, &mut tq), "run quality")?;
            }
            record(self, t as i64)?;
            unsafe {
                check(sys::EN_nextH(self.ph, &mut step), "advance hydraulics")?;
                check(sys::EN_nextQ(self.ph, &mut qstep), "advance quality")?;
            }
            if step <= 0 {
                return Ok(());
            }
        })();
        unsafe {
            sys::EN_closeQ(self.ph);
            sys::EN_closeH(self.ph);
        }
        result
    }
}

impl Drop for Project {
    fn drop(&mut self) {
        // SAFETY: ph came from EN_createproject and is not used after this.
        unsafe {
            sys::EN_close(self.ph);
            sys::EN_deleteproject(self.ph);
        }
        for f in &self.temp_files {
            let _ = std::fs::remove_file(f);
        }
    }
}
