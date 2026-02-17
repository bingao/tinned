use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

use tinned::TinnedError;

// Python exception type
create_exception!(tinned, TinnedErrorPy, PyException);

// Local wrapper type
pub struct PyTinnedError(TinnedError);

impl From<TinnedError> for PyTinnedError {
    fn from(err: TinnedError) -> Self {
        Self(err)
    }
}

impl From<PyTinnedError> for PyErr {
    fn from(err: PyTinnedError) -> Self {
        let msg = err.0.to_string();
        // Create the Python exception with (msg,) as args.
        PyErr::new::<TinnedErrorPy, _>((msg,))
    }
}

// Convenience helper so you can write: return Err(to_pyerr(err));
pub fn to_pyerr(err: TinnedError) -> PyErr {
    PyTinnedError::from(err).into()
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("TinnedError", m.py().get_type::<TinnedErrorPy>())?;
    Ok(())
}
