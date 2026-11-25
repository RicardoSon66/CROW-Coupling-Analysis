"""
2-ring CROW (add-drop) 해석 모델을 FDTD drop 스펙트럼에 피팅
  - 입력: data/crow_gap*.csv  (export_crow_csv.py 가 만든 파일)
  - 출력: data/tmm_fit_params.csv, data/tmm_fit.png, data/tmm_gap_design.png

모델 (Mason 신호흐름 그래프로 유도, 동일한 두 링):
  t_i^2 + k_i^2 = 1,  r = a * exp(-j*phi)  (링 1바퀴),  h = sqrt(r) (반 바퀴)
  k1 = k3 = k_bus (링-버스),  k2 = k_rr (링-링)

  Delta = 1 - t2 (t1 + t3) r + t1 t3 r^2
  Drop  = -j k1 k2 k3 r / Delta
  Thru  = [ t1 Delta - k1^2 t2 r (1 - t2 t3 r) + k1^2 k2^2 t3 r^2 ] / Delta

  phi(lam) = 2*pi * n_g * L * (1/lam - 1/lam0)   (lam0 = 단일 링 공진 파장)

피팅 파라미터 (세 gap 동시 피팅):
  공통: lam0, a(링 1바퀴 진폭 투과), k_bus
  gap별: k_rr
  n_g 는 넓은 대역 CSV의 FSR에서 구해 고정 (좁은 창에서는 k 와 n_g 가 서로 구분 안 됨)
"""
import glob
import re
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_UM = 4.0                       # 링 중심선 반지름 (um)
L_UM = 2 * np.pi * R_UM          # 링 둘레
LAM_MIN, LAM_MAX = 1.546, 1.552  # 피팅 창
FLOOR_DB = -30.0                 # 이보다 낮은 값은 잘라서 비교 (truncation 물결 무시)
NG_FALLBACK = 4.2


# ----------------------------------------------------------------- 모델
def crow2(lam, lam0, ng, a, k_bus, k_rr):
    phi = 2 * np.pi * ng * L_UM * (1.0 / lam - 1.0 / lam0)
    r = a * np.exp(-1j * phi)
    k1 = k3 = k_bus
    k2 = k_rr
    t1 = t3 = np.sqrt(1 - k_bus ** 2)
    t2 = np.sqrt(1 - k_rr ** 2)
    delta = 1 - t2 * (t1 + t3) * r + t1 * t3 * r ** 2
    drop = -1j * k1 * k2 * k3 * r / delta
    thru = (t1 * delta - k1 ** 2 * t2 * r * (1 - t2 * t3 * r)
            + k1 ** 2 * k2 ** 2 * t3 * r ** 2) / delta
    return np.abs(thru) ** 2, np.abs(drop) ** 2


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-12))


# ----------------------------------------------------------- 데이터 읽기
def load_csvs():
    sets = []
    for f in sorted(glob.glob("data/crow_gap*.csv")):
        m = re.search(r"gap([0-9.]+?)(?:_|\.csv)", f)
        gap = float(m.group(1))
        df = pd.read_csv(f)
        sets.append(dict(file=f, gap=gap, lam=df["wavelength_um"].values,
                         drop_db=df["drop_dB"].values))
    if not sets:
        raise SystemExit("data/crow_gap*.csv 가 없습니다. export_crow_csv.py 를 먼저 실행하세요.")
    return sets


def ng_from_fsr(sets):
    """넓은 대역(>30 nm) CSV에서 공진 중심들을 찾아 FSR -> n_g"""
    for s in sets:
        lam, d = s["lam"], s["drop_db"]
        if lam.max() - lam.min() < 0.030:
            continue
        on = d > -10
        groups, cur = [], []
        for i, flag in enumerate(on):
            if flag:
                cur.append(i)
            elif cur:
                groups.append(cur)
                cur = []
        if cur:
            groups.append(cur)
        centers = [np.average(lam[g], weights=10 ** (d[g] / 10)) for g in groups]
        if len(centers) >= 2:
            centers = np.array(centers)
            fsr = np.diff(centers)
            mid = 0.5 * (centers[1:] + centers[:-1])
            i = np.argmin(np.abs(mid - 1.5492))
            ng = mid[i] ** 2 / (fsr[i] * L_UM)
            print(f"n_g from FSR ({s['file']}): FSR={fsr[i]*1e3:.2f} nm -> n_g={ng:.4f}")
            return ng
    print(f"넓은 대역 CSV가 없어 n_g={NG_FALLBACK} 사용")
    return NG_FALLBACK


# ------------------------------------------------------------------ 피팅
def fit(sets, ng):
    win = []
    for s in sets:
        m = (s["lam"] >= LAM_MIN) & (s["lam"] <= LAM_MAX)
        win.append((s["lam"][m], np.maximum(s["drop_db"][m], FLOOR_DB)))

    n = len(sets)
    # p = [lam0, a, k_bus, k_rr_1..k_rr_n]
    p0 = [1.54921, 0.995, 0.25] + [0.08 * np.exp(-(s["gap"] - 0.35) / 0.07) for s in sets]
    lo = [1.5487, 0.90, 0.02] + [1e-4] * n
    hi = [1.5497, 1.00, 0.80] + [0.80] * n

    def resid(p):
        out = []
        for i, (lam, d) in enumerate(win):
            _, D = crow2(lam, p[0], ng, p[1], p[2], p[3 + i])
            out.append(np.maximum(db(D), FLOOR_DB) - d)
        return np.concatenate(out)

    res = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac")
    rms = np.sqrt(np.mean(res.fun ** 2))
    return res.x, rms


# --------------------------------------------------- 통과대역 지표 (모델)
def model_metrics(lam0, ng, a, k_bus, k_rr):
    lam = np.linspace(lam0 - 0.0015, lam0 + 0.0015, 30001)
    _, D = crow2(lam, lam0, ng, a, k_bus, k_rr)
    d = db(D)
    idx = np.where((d[1:-1] > d[:-2]) & (d[1:-1] > d[2:]))[0] + 1
    top = np.sort(idx[np.argsort(d[idx])[-2:]]) if len(idx) >= 2 else idx
    ripple = 0.0
    if len(top) == 2:
        dip = d[top[0]:top[1] + 1].min()
        ripple = max(d[top].mean() - dip, 0.0)
    ref = d.max()
    above = np.where(d >= ref - 3)[0]
    bw = (lam[above[-1]] - lam[above[0]]) * 1e3
    return ripple, bw, -ref


def k_flat(lam0, ng, a, k_bus):
    """리플이 사라지는 경계의 k_rr (maximally flat 근사)"""
    f = lambda k: model_metrics(lam0, ng, a, k_bus, k)[0] - 0.01
    ks = np.geomspace(1e-3, 0.3, 80)          # 약한 결합 -> 강한 결합 순으로 스캔
    vals = [f(k) for k in ks]
    for i in range(1, len(ks)):
        if vals[i - 1] < 0 <= vals[i]:        # 리플이 처음 생기는 경계
            return brentq(f, ks[i - 1], ks[i], xtol=1e-6)
    raise RuntimeError("스캔 범위에서 리플 경계를 찾지 못했습니다.")


# ------------------------------------------------------------------ main
def main():
    sets = load_csvs()
    ng = ng_from_fsr(sets)
    p, rms = fit(sets, ng)
    lam0, a, k_bus = p[:3]
    k_rr = p[3:]
    gaps = np.array([s["gap"] for s in sets])

    loss_db_rt = -20 * np.log10(a)
    alpha_db_cm = loss_db_rt / (L_UM * 1e-4)
    print("\n=== 피팅 결과 (RMS 오차 %.2f dB) ===" % rms)
    print(f"lam0 = {lam0:.6f} um,  n_g = {ng:.4f}")
    print(f"a = {a:.5f}  (1바퀴 손실 {loss_db_rt:.3f} dB, 약 {alpha_db_cm:.1f} dB/cm)")
    print(f"k_bus = {k_bus:.4f}  (power {k_bus**2:.4f})")

    rows = []
    for g, k in zip(gaps, k_rr):
        rip, bw, il = model_metrics(lam0, ng, a, k_bus, k)
        rows.append(dict(gap_um=g, k_rr=k, k_rr_power=k ** 2,
                         model_ripple_dB=rip, model_bw3dB_nm=bw, model_IL_dB=il))
        print(f"gap {g:.2f} um: k_rr = {k:.4f} (power {k**2:.5f}), "
              f"ripple {rip:.2f} dB, BW {bw:.3f} nm, IL {il:.2f} dB")

    # k_rr(gap) = k0 * exp(-gamma * gap)  (지수 감쇠 피팅)
    c1, c0 = np.polyfit(gaps, np.log(k_rr), 1)
    gamma = -c1
    print(f"\nk_rr(gap) ~ {np.exp(c0):.3f} * exp(-{gamma:.2f} * gap[um])"
          f"  (0.1 um 당 field x{np.exp(-gamma*0.1):.3f})")

    kf = k_flat(lam0, ng, a, k_bus)
    g_flat = (np.log(kf) - c0) / c1
    rip_f, bw_f, il_f = model_metrics(lam0, ng, a, k_bus, kf)
    print(f"\n[예측] 리플이 사라지는 k_rr = {kf:.4f} (power {kf**2:.5f})")
    print(f"[예측] 해당 gap = {g_flat:.3f} um, BW3dB = {bw_f:.3f} nm "
          f"({bw_f*124.9:.1f} GHz), IL = {il_f:.2f} dB")
    print(f"참고: 근사식 k_rr^2 = (k_bus^2)^2/4 -> {k_bus**4/4:.5f}")

    out = pd.DataFrame(rows)
    out["lam0_um"], out["n_g"], out["a"], out["k_bus"] = lam0, ng, a, k_bus
    out["pred_flat_gap_um"], out["pred_flat_k_rr"] = g_flat, kf
    out.to_csv("data/tmm_fit_params.csv", index=False)

    # ---- 그림 1: FDTD vs TMM
    fig, axes = plt.subplots(1, len(sets), figsize=(5 * len(sets), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, s, k in zip(axes, sets, k_rr):
        m = (s["lam"] >= 1.5482) & (s["lam"] <= 1.5502)
        lam_f = np.linspace(1.5482, 1.5502, 4001)
        _, D = crow2(lam_f, lam0, ng, a, k_bus, k)
        ax.plot((s["lam"][m] - lam0) * 1e3, s["drop_db"][m], "o", ms=2.5, label="FDTD")
        ax.plot((lam_f - lam0) * 1e3, db(D), "-", lw=1.5, label="TMM fit")
        ax.set_title(f"gap = {s['gap']:.2f} µm  (κ_rr = {k:.3f})")
        ax.set_xlabel(f"λ − {lam0*1e3:.3f} nm  (nm)")
        ax.set_ylim(-30, 1)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Drop (dB)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig("data/tmm_fit.png", dpi=200)

    # ---- 그림 2: gap 설계 곡선 (리플 / 대역폭)
    g_axis = np.linspace(gaps.min() - 0.05, max(gaps.max(), g_flat) + 0.1, 61)
    rip_c, bw_c = [], []
    for g in g_axis:
        k = np.exp(c0 + c1 * g)
        r_, b_, _ = model_metrics(lam0, ng, a, k_bus, k)
        rip_c.append(r_)
        bw_c.append(b_)
    fdtd = pd.read_csv("data/crow_summary.csv") if glob.glob("data/crow_summary.csv") else None
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(g_axis, rip_c, "-", color="C0", label="ripple (TMM)")
    if fdtd is not None:
        ax1.plot(fdtd["gap_um"], fdtd["ripple_dB"], "o", color="C0", label="ripple (FDTD)")
    ax1.axvline(g_flat, ls="--", color="gray")
    ax1.set_xlabel("Ring-ring gap (µm)")
    ax1.set_ylabel("Ripple (dB)", color="C0")
    ax2 = ax1.twinx()
    ax2.plot(g_axis, np.array(bw_c) * 124.9, "-", color="C1", label="BW (TMM)")
    if fdtd is not None:
        ax2.plot(fdtd["gap_um"], fdtd["bw3dB_nm"] * 124.9, "s", color="C1", label="BW (FDTD)")
    ax2.set_ylabel("3 dB bandwidth (GHz)", color="C1")
    ax1.set_title(f"Predicted flat-top gap ≈ {g_flat:.3f} µm")
    fig.legend(loc="upper right", bbox_to_anchor=(0.88, 0.88), fontsize=8)
    fig.tight_layout()
    fig.savefig("data/tmm_gap_design.png", dpi=200)
    print("\n저장: data/tmm_fit_params.csv, data/tmm_fit.png, data/tmm_gap_design.png")


if __name__ == "__main__":
    main()
