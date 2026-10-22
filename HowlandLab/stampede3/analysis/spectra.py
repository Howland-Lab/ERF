#!/usr/bin/env python3
"""spectra.py -- the lab's spectrum conventions, written down once (22 Oct 2026) so that every spectrum we show is
computed by a function that has passed the known-answer tests in verify_analysis.py.

The 15 Oct tile-imprint figure (wrf_tile_imprint_2panel.png) and the tile-mode shares (scratchpad tilespec.py) used the
2-D convention below: |FFT2(field - mean)|^2 with the DC mode zeroed, summed over frames.

Conventions
  spectrum1d(x, dx) -> k, E     one-sided, k = 2 pi m / L (rad/m), m = 1 .. n//2, DC removed, normalised so that
                                sum(E) * dk = var(x)  (Parseval), dk = 2 pi / L.  E has units of x^2 per (rad/m).
  spectrum2d(a) -> P            P[ky, kx] = |FFT2(a - mean)|^2 / N^2 in numpy FFT ordering (use np.fft.fftfreq(n) * n for
                                the integer wavenumbers), P[0, 0] = 0, sum(P) = var(a)  (Parseval).
  line_share(P, py, axis=0)     fraction of the 2-D energy on the rows ky = n * py (n != 0): a field made of ny/py
                                identical copies along y has all its energy there (share = 1); a field without that
                                periodicity has share ~ 1/py.  This is the WRF tiling diagnosis.
  tile_mode_share(P, px, py)    fraction of the energy in the modes (m px, n py), |m|,|n| <= 3, as in tilespec.py.
"""
import numpy as np


def spectrum1d(x, dx):
    x = np.asarray(x, float); x = x - x.mean(); n = x.size; L = n * dx
    X = np.fft.rfft(x)
    P = np.abs(X) ** 2 / n ** 2                 # two-sided power per mode; sum over all n modes = var(x)
    E = 2.0 * P[1:]                             # fold negative frequencies onto m = 1 .. n//2
    if n % 2 == 0:
        E[-1] = P[-1]                           # the Nyquist mode appears once
    m = np.arange(1, E.size + 1); dk = 2 * np.pi / L
    return m * dk, E / dk


def spectrum2d(a):
    a = np.asarray(a, float); a = a - a.mean()
    P = np.abs(np.fft.fft2(a)) ** 2 / a.size ** 2
    P[0, 0] = 0.0
    return P


def line_share(P, py, axis=0, include_zero=False):
    """Energy on the lines k = n * py along `axis`, as a fraction of the total.  By default the n = 0 line is left out
    (that is the tiling diagnosis: a field made of ny/py copies has share (py-1)/py there, an ordinary field ~ 1/py);
    with include_zero=True a perfectly periodic field gives exactly 1."""
    nk = P.shape[axis]; k = (np.fft.fftfreq(nk) * nk).round().astype(int)
    on = (k % py == 0) & ((k != 0) | include_zero)
    sel = np.take(P, np.where(on)[0], axis=axis)
    return float(sel.sum() / P.sum())


def tile_mode_share(P, px, py, mmax=3):
    ny, nx = P.shape
    tile = sum(P[(n * py) % ny, (m * px) % nx] for m in range(-mmax, mmax + 1) for n in range(-mmax, mmax + 1) if (m, n) != (0, 0))
    return float(tile / P.sum())
