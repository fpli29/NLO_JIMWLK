"""Non-production coordinate-space KLM kernel diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .coordinate_kernels import dot, norm2, safe_inv, validate_coords, vec


IMPLEMENTED_UNBARRED_KERNELS = ("KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ")
PENDING_UNBARRED_KERNELS = ()


@dataclass(frozen=True)
class KJSJIntegrationPolicy:
    """Explicit finite-grid policy for the K_JSJ tilde-K integral.

    This policy is diagnostic bookkeeping, not a physical regulator.
    """

    quadrature_weights: np.ndarray
    mu: float
    excluded_indices: tuple[int, ...] = ()
    exclude_coincident_labels: tuple[str, ...] = ()
    principal_value: str = "none"
    subtraction: str = "diagonal_zero"
    finite_volume_boundary: str = "finite_coordinate_sum"
    description: str = "explicit diagnostic finite z-prime sum"

    def validate(self, ncoords: int) -> np.ndarray:
        weights = np.asarray(self.quadrature_weights, dtype=float)
        if weights.shape != (ncoords,):
            raise ValueError(f"quadrature_weights must have shape {(ncoords,)}, got {weights.shape}")
        if not np.all(np.isfinite(weights)):
            raise ValueError("quadrature_weights must be finite")
        if self.mu <= 0.0 or not np.isfinite(self.mu):
            raise ValueError("K_JSJ policy mu must be finite and positive")
        allowed_labels = {"x", "y", "z"}
        unknown_labels = set(self.exclude_coincident_labels) - allowed_labels
        if unknown_labels:
            raise ValueError(f"unknown exclude_coincident_labels: {sorted(unknown_labels)}")
        if self.principal_value != "none":
            raise NotImplementedError("only principal_value='none' is implemented for diagnostics")
        if self.subtraction not in {"none", "diagonal_zero"}:
            raise NotImplementedError("only subtraction='none' or 'diagonal_zero' is implemented")
        if self.finite_volume_boundary not in {"finite_coordinate_sum", "periodic_box_label_only"}:
            raise NotImplementedError(
                "finite_volume_boundary must be 'finite_coordinate_sum' or 'periodic_box_label_only'"
            )
        for index in self.excluded_indices:
            if index < 0 or index >= ncoords:
                raise ValueError(f"excluded index {index} is out of range for {ncoords} coordinates")
        return weights

    def excluded_for(self, x: int, y: int, z: int, ncoords: int) -> set[int]:
        self.validate(ncoords)
        excluded = set(int(i) for i in self.excluded_indices)
        label_to_index = {"x": int(x), "y": int(y), "z": int(z)}
        for label in self.exclude_coincident_labels:
            excluded.add(label_to_index[label])
        return excluded


def _as_complex_if_needed(value):
    return np.real_if_close(value)


def _safe_positive(value: float, *, singularity_policy: str, eps: float | None, name: str) -> float:
    val = float(value)
    if val > 0.0:
        return val
    if singularity_policy == "raise":
        raise ValueError(f"singular {name}: logarithm argument is nonpositive")
    if singularity_policy == "nan":
        return float("nan")
    if singularity_policy == "eps":
        if eps is None or eps <= 0.0:
            raise ValueError("eps policy requires a positive eps argument")
        return float(eps)
    raise ValueError("singularity_policy must be 'raise', 'nan', or 'eps'")


def _log_ratio(num: float, den: float, *, singularity_policy: str, eps: float | None, name: str) -> float:
    safe_num = _safe_positive(num, singularity_policy=singularity_policy, eps=eps, name=f"{name} numerator")
    safe_den = _safe_positive(den, singularity_policy=singularity_policy, eps=eps, name=f"{name} denominator")
    return float(np.log(safe_num / safe_den))


def _inv(value: float, *, singularity_policy: str, eps: float | None, name: str) -> float:
    return safe_inv(value, singularity_policy=singularity_policy, eps=eps, name=name)


def _geom_x_y_z_zp(coords, x, y, z, zp):
    X = vec(coords, x, z)
    Xp = vec(coords, x, zp)
    Y = vec(coords, y, z)
    Yp = vec(coords, y, zp)
    rxy = vec(coords, x, y)
    rzzp = vec(coords, z, zp)
    return {
        "X": X,
        "Xp": Xp,
        "Y": Y,
        "Yp": Yp,
        "rxy": rxy,
        "rzzp": rzzp,
        "X2": norm2(X),
        "Xp2": norm2(Xp),
        "Y2": norm2(Y),
        "Yp2": norm2(Yp),
        "rxy2": norm2(rxy),
        "rzzp2": norm2(rzzp),
    }


def KJSJ_unbarred_value(
    coords,
    x,
    y,
    z,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
    integration_policy: KJSJIntegrationPolicy | None = None,
):
    """Return diagnostic K_JSJ(x,y;z) from WORKNLO.tex lines 324--332."""

    if integration_policy is None:
        raise ValueError("K_JSJ requires an explicit KJSJIntegrationPolicy")
    if int(x) == int(y) and integration_policy.subtraction == "diagonal_zero":
        return 0.0
    local = KJSJ_unbarred_local_value(
        coords,
        x,
        y,
        z,
        Nc=Nc,
        nf=nf,
        alpha_s=alpha_s,
        mu=integration_policy.mu,
        singularity_policy=singularity_policy,
        eps=eps,
    )
    tilde_integral = KJSJ_unbarred_tilde_integral_value(
        coords,
        x,
        y,
        z,
        Nc=Nc,
        nf=nf,
        alpha_s=alpha_s,
        singularity_policy=singularity_policy,
        eps=eps,
        integration_policy=integration_policy,
    )
    return _as_complex_if_needed(local - 0.5 * Nc * tilde_integral)


def KJSJ_beta0(Nc=3, nf=0) -> float:
    """Return b = 11/3 Nc - 2/3 nf from WORKNLO.tex line 332."""

    return (11.0 / 3.0) * Nc - (2.0 / 3.0) * nf


def KJSJ_unbarred_local_value(
    coords,
    x,
    y,
    z,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    mu=1.0,
    singularity_policy="raise",
    eps=None,
):
    """Return the local coordinate/scheme part of K_JSJ before the tilde integral."""

    coords = validate_coords(coords)
    if int(x) == int(y):
        return 0.0
    if mu <= 0.0 or not np.isfinite(mu):
        raise ValueError("mu must be finite and positive")
    X = vec(coords, x, z)
    Y = vec(coords, y, z)
    rxy = vec(coords, x, y)
    X2 = norm2(X)
    Y2 = norm2(Y)
    rxy2 = norm2(rxy)
    inv_X2 = _inv(X2, singularity_policy=singularity_policy, eps=eps, name="X^2")
    inv_Y2 = _inv(Y2, singularity_policy=singularity_policy, eps=eps, name="Y^2")
    inv_rxy2 = _inv(rxy2, singularity_policy=singularity_policy, eps=eps, name="(x-y)^2")
    b0 = KJSJ_beta0(Nc=Nc, nf=nf)
    scheme_log = _log_ratio(
        rxy2 * mu**2,
        1.0,
        singularity_policy=singularity_policy,
        eps=eps,
        name="(x-y)^2 mu^2",
    )
    xy_log = _log_ratio(X2, Y2, singularity_policy=singularity_policy, eps=eps, name="X^2/Y^2")
    constant = ((67.0 / 9.0) - (np.pi**2 / 3.0)) * Nc - (10.0 / 9.0) * nf
    bracket = b0 * scheme_log - b0 * (X2 - Y2) * inv_rxy2 * xy_log + constant
    return -alpha_s**2 / (16.0 * np.pi**3) * rxy2 * inv_X2 * inv_Y2 * bracket


def KJSJ_unbarred_tilde_integral_value(
    coords,
    x,
    y,
    z,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
    integration_policy: KJSJIntegrationPolicy,
):
    """Return the explicit finite-grid diagnostic sum for int dz' tilde_K."""

    coords = validate_coords(coords)
    ncoords = coords.shape[0]
    weights = integration_policy.validate(ncoords)
    if int(x) == int(y) and integration_policy.subtraction == "diagonal_zero":
        return 0.0
    excluded = integration_policy.excluded_for(x, y, z, ncoords)
    total = 0.0 + 0.0j
    for zp, weight in enumerate(weights):
        if zp in excluded or weight == 0.0:
            continue
        total += weight * tilde_K_JJSSJ_unbarred_value(
            coords,
            x,
            y,
            z,
            zp,
            Nc=Nc,
            nf=nf,
            alpha_s=alpha_s,
            singularity_policy=singularity_policy,
            eps=eps,
        )
    return _as_complex_if_needed(total)


def KJJSSJ_unbarred_value(
    coords,
    w,
    x,
    y,
    z,
    zp,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
):
    """Return K_JJSSJ(w;x,y;z,z') from WORKNLO.tex lines 288--297."""

    _ = (Nc, nf)
    coords = validate_coords(coords)
    g = _geom_x_y_z_zp(coords, x, y, z, zp)
    W = vec(coords, w, z)
    Wp = vec(coords, w, zp)
    zpmz = vec(coords, zp, z)
    zmzp = vec(coords, z, zp)
    W2 = norm2(W)
    Wp2 = norm2(Wp)
    inv_X2 = _inv(g["X2"], singularity_policy=singularity_policy, eps=eps, name="X^2")
    inv_Yp2 = _inv(g["Yp2"], singularity_policy=singularity_policy, eps=eps, name="Y'^2")
    inv_W2 = _inv(W2, singularity_policy=singularity_policy, eps=eps, name="W^2")
    inv_Wp2 = _inv(Wp2, singularity_policy=singularity_policy, eps=eps, name="W'^2")
    inv_zzp2 = _inv(g["rzzp2"], singularity_policy=singularity_policy, eps=eps, name="(z-z')^2")

    tensor = np.outer(g["X"], g["Yp"]) * inv_X2 * inv_Yp2
    bracket = (
        np.eye(2) * (0.5 * inv_zzp2)
        + np.outer(zpmz, Wp) * inv_zzp2 * inv_Wp2
        + np.outer(W, zmzp) * inv_zzp2 * inv_W2
        - np.outer(W, Wp) * inv_W2 * inv_Wp2
    )
    log_term = _log_ratio(
        W2,
        Wp2,
        singularity_policy=singularity_policy,
        eps=eps,
        name="W^2/W'^2",
    )
    contraction = float(np.sum(tensor * bracket))
    return -1.0j * alpha_s**2 / (2.0 * np.pi**4) * contraction * log_term


def tilde_K_JJSSJ_unbarred_value(
    coords,
    x,
    y,
    z,
    zp,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
):
    """Return tilde K from WORKNLO.tex lines 307--311."""

    params = {
        "Nc": Nc,
        "nf": nf,
        "alpha_s": alpha_s,
        "singularity_policy": singularity_policy,
        "eps": eps,
    }
    return 0.5j * (
        KJJSSJ_unbarred_value(coords, x, x, y, z, zp, **params)
        - KJJSSJ_unbarred_value(coords, y, x, y, z, zp, **params)
        - KJJSSJ_unbarred_value(coords, x, y, x, z, zp, **params)
        + KJJSSJ_unbarred_value(coords, y, y, x, z, zp, **params)
    )


def KJJSJ_unbarred_value(
    coords,
    w,
    x,
    y,
    z,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
):
    """Return K_JJSJ(w;x,y;z) from WORKNLO.tex lines 298--300."""

    _ = (Nc, nf)
    coords = validate_coords(coords)
    if int(x) == int(y):
        return 0.0
    X = vec(coords, x, z)
    Y = vec(coords, y, z)
    W = vec(coords, w, z)
    rxy = vec(coords, x, y)
    X2 = norm2(X)
    Y2 = norm2(Y)
    W2 = norm2(W)
    rxy2 = norm2(rxy)
    term = (
        dot(X, W)
        * _inv(X2, singularity_policy=singularity_policy, eps=eps, name="X^2")
        * _inv(W2, singularity_policy=singularity_policy, eps=eps, name="W^2")
        - dot(Y, W)
        * _inv(Y2, singularity_policy=singularity_policy, eps=eps, name="Y^2")
        * _inv(W2, singularity_policy=singularity_policy, eps=eps, name="W^2")
    )
    log_y = _log_ratio(Y2, rxy2, singularity_policy=singularity_policy, eps=eps, name="Y^2/(x-y)^2")
    log_x = _log_ratio(X2, rxy2, singularity_policy=singularity_policy, eps=eps, name="X^2/(x-y)^2")
    return -1.0j * alpha_s**2 / (4.0 * np.pi**3) * term * log_y * log_x


def Kqbarq_unbarred_value(
    coords,
    x,
    y,
    z,
    zp,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
):
    """Return K_qbarq(x,y;z,z') from WORKNLO.tex lines 301--306."""

    _ = Nc
    coords = validate_coords(coords)
    if int(x) == int(y):
        return 0.0
    g = _geom_x_y_z_zp(coords, x, y, z, zp)
    delta = g["X2"] * g["Yp2"] - g["Xp2"] * g["Y2"]
    numerator = g["Xp2"] * g["Y2"] + g["Yp2"] * g["X2"] - g["rxy2"] * g["rzzp2"]
    inv_zzp4 = _inv(g["rzzp2"] ** 2, singularity_policy=singularity_policy, eps=eps, name="(z-z')^4")
    inv_delta = _inv(delta, singularity_policy=singularity_policy, eps=eps, name="X^2Y'^2-X'^2Y^2")
    log_term = _log_ratio(
        g["X2"] * g["Yp2"],
        g["Xp2"] * g["Y2"],
        singularity_policy=singularity_policy,
        eps=eps,
        name="X^2Y'^2/(X'^2Y^2)",
    )
    bracket = numerator * inv_zzp4 * inv_delta * log_term - 2.0 * inv_zzp4
    return -alpha_s**2 * nf / (8.0 * np.pi**4) * bracket


def KJSSJ_unbarred_value(
    coords,
    x,
    y,
    z,
    zp,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
):
    """Return K_JSSJ(x,y;z,z') from WORKNLO.tex lines 313--323."""

    coords = validate_coords(coords)
    if int(x) == int(y):
        return 0.0
    g = _geom_x_y_z_zp(coords, x, y, z, zp)
    delta = g["X2"] * g["Yp2"] - g["Xp2"] * g["Y2"]
    inv_zzp2 = _inv(g["rzzp2"], singularity_policy=singularity_policy, eps=eps, name="(z-z')^2")
    inv_zzp4 = _inv(g["rzzp2"] ** 2, singularity_policy=singularity_policy, eps=eps, name="(z-z')^4")
    inv_delta = _inv(delta, singularity_policy=singularity_policy, eps=eps, name="X^2Y'^2-X'^2Y^2")
    inv_X2Yp2 = _inv(g["X2"] * g["Yp2"], singularity_policy=singularity_policy, eps=eps, name="X^2Y'^2")
    inv_Y2Xp2 = _inv(g["Y2"] * g["Xp2"], singularity_policy=singularity_policy, eps=eps, name="Y^2X'^2")
    inv_Xp2Y2 = _inv(g["Xp2"] * g["Y2"], singularity_policy=singularity_policy, eps=eps, name="X'^2Y^2")
    log_term = _log_ratio(
        g["X2"] * g["Yp2"],
        g["Xp2"] * g["Y2"],
        singularity_policy=singularity_policy,
        eps=eps,
        name="X^2Y'^2/(X'^2Y^2)",
    )
    curly = (
        2.0
        * (g["X2"] * g["Yp2"] + g["Xp2"] * g["Y2"] - 4.0 * g["rxy2"] * g["rzzp2"])
        * inv_zzp4
        * inv_delta
        + g["rxy2"] ** 2 * inv_delta * (inv_X2Yp2 + inv_Y2Xp2)
        + g["rxy2"] * inv_zzp2 * (inv_X2Yp2 - inv_Xp2Y2)
    )
    base = alpha_s**2 / (16.0 * np.pi**4) * (-4.0 * inv_zzp4 + curly * log_term)
    tilde = tilde_K_JJSSJ_unbarred_value(
        coords,
        x,
        y,
        z,
        zp,
        Nc=Nc,
        nf=nf,
        alpha_s=alpha_s,
        singularity_policy=singularity_policy,
        eps=eps,
    )
    return _as_complex_if_needed(base + tilde)


def _build_array(coords, shape, value_func, *, singularity_policy="raise", eps=None, **params):
    validate_coords(coords)
    out = np.empty(shape, dtype=complex)
    for index in np.ndindex(shape):
        out[index] = value_func(
            coords,
            *index,
            singularity_policy=singularity_policy,
            eps=eps,
            **params,
        )
    return np.real_if_close(out)


def build_KJSJ_unbarred(coords, **params):
    """Return array shape (N,N,N): K[x,y,z]."""

    n = validate_coords(coords).shape[0]
    return _build_array(coords, (n, n, n), KJSJ_unbarred_value, **params)


def build_KJSSJ_unbarred(coords, **params):
    """Return array shape (N,N,N,N): K[x,y,z,zp]."""

    n = validate_coords(coords).shape[0]
    return _build_array(coords, (n, n, n, n), KJSSJ_unbarred_value, **params)


def build_Kqbarq_unbarred(coords, **params):
    """Return array shape (N,N,N,N): K[x,y,z,zp]."""

    n = validate_coords(coords).shape[0]
    return _build_array(coords, (n, n, n, n), Kqbarq_unbarred_value, **params)


def build_KJJSJ_unbarred(coords, **params):
    """Return array shape (N,N,N,N): K[w,x,y,z]."""

    n = validate_coords(coords).shape[0]
    return _build_array(coords, (n, n, n, n), KJJSJ_unbarred_value, **params)


def build_KJJSSJ_unbarred(coords, **params):
    """Return array shape (N,N,N,N,N): K[w,x,y,z,zp]."""

    n = validate_coords(coords).shape[0]
    return _build_array(coords, (n, n, n, n, n), KJJSSJ_unbarred_value, **params)


def build_all_unbarred_physical_kernels(coords, **params):
    """Build implemented unbarred physical kernels plus metadata."""

    kjsj_params = dict(params)
    shared_params = dict(params)
    shared_params.pop("integration_policy", None)
    kernels = {
        "KJSJ": build_KJSJ_unbarred(coords, **kjsj_params),
        "KJSSJ": build_KJSSJ_unbarred(coords, **shared_params),
        "Kqbarq": build_Kqbarq_unbarred(coords, **shared_params),
        "KJJSJ": build_KJJSJ_unbarred(coords, **shared_params),
        "KJJSSJ": build_KJJSSJ_unbarred(coords, **shared_params),
    }
    kernels["metadata"] = {
        "implemented_kernels": list(IMPLEMENTED_UNBARRED_KERNELS),
        "pending_kernels": list(PENDING_UNBARRED_KERNELS),
        "nonproduction_only": True,
        "kernel_type": "unbarred singlet",
    }
    return kernels
