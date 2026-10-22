#!/usr/bin/env python3
"""planes_clean_fig.py -- deck version of planes_w_cbl_sd25.png without the invalid tiled WRF run (from planes.json)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, "figs")
b = np.linspace(-8, 8, 321); c = 0.5 * (b[1:] + b[:-1]); dw = b[1] - b[0]          # HBINS["w"] in erf_planes.py
R = 4                                                   # merge 4 native bins (0.2 m/s): fewer single-count bins in the tails
cr = c.reshape(-1, R).mean(1); dwr = dw * R


def merged_pdf(hist, n, std, bins=b, R=R, min_count=5):
    """Merge R native histogram bins and hide bins with fewer than min_count samples (sampling noise).
    Returns (w'/std at the merged-bin centres, pdf of w'/std as a masked array).  Normalisation: sum(pdf) * (R dw / std)
    = fraction of the n samples that fall inside the bins (= 1 when none are outside).  Factored out 22 Oct 2026 for
    verify_analysis.py; identical to the inline code used for the deck figure."""
    c = 0.5 * (bins[1:] + bins[:-1]); dw = bins[1] - bins[0]
    cr = c.reshape(-1, R).mean(1); dwr = dw * R
    h = np.asarray(hist, float).reshape(-1, R).sum(1)
    pdf = h / n / dwr * std
    return cr / std, np.ma.masked_where(h < min_count, pdf)


def curve(q):
    return merged_pdf(q["hist"], q["n"], q["std"])


def main():
  D = json.load(open(os.path.join(HERE, "planes.json")))["cbl_sd25"]
  E, W = D["ERF"], D["WRF"]
  labs = [k[:-2] for k in E if k.endswith("_w")]
  fig, ax = plt.subplots(1, len(labs), figsize=(3.6 * len(labs), 3.4), sharey=True)
  for j, lab in enumerate(labs):
      e = E[f"{lab}_w"]
      x, y = curve(e); ax[j].semilogy(x, y, color="#9B2335", lw=1.6, label=f"ERF 25 m (N = {e['n'] / 1e6:.1f} M)")
      for w, col in (("cbl_dx30", "#555555"), ("cbl_dx10", "#222222")):
          q = W.get(w, {}).get(f"{lab}_w")
          if q is not None:
              x, y = curve(q); ax[j].semilogy(x, y, color=col, lw=1.0, label=f"WRF {w.split('_')[1][2:]} m (N = {q['n'] / 1e6:.1f} M)")
      ax[j].semilogy(c, np.exp(-0.5 * c ** 2) / np.sqrt(2 * np.pi), "--", color="#999999", lw=0.8, label="Gaussian")
      ax[j].set_xlim(-5, 7); ax[j].set_ylim(1e-6, 1); ax[j].grid(alpha=0.25)
      ax[j].set_title(f"z/zi = {lab}" if lab != "hub" else "z = 100 m", fontsize=10); ax[j].set_xlabel(r"$w'/\sigma_w$")
  ax[0].set_ylabel("PDF"); ax[0].legend(frameon=False, fontsize=6.5, loc="lower center")
  fig.suptitle("Vertical-velocity distributions, convective case, last hour", fontsize=11)
  fig.tight_layout(); fig.savefig(f"{FIG}/planes_w_cbl_sd25_clean.png", dpi=150); plt.close(fig)
  print("labs", labs)


if __name__ == "__main__":
    main()
