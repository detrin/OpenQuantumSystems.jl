
const ComputableType = Union{AbstractFloat,Complex}
const version = "0.2.0"

const _SAFE_DIV_TOL = 1e-10

_safe_div(a, b; tol=_SAFE_DIV_TOL) = abs(b) < tol ? zero(a) : a / b

_safe_inv(x; tol=_SAFE_DIV_TOL) = abs(x) < tol ? zero(x) : 1 / x