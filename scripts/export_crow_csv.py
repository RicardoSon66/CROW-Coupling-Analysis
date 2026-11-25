"""
CROW gap sweep: 같은 폴더의 Tidy3D 결과 hdf5 -> CSV 저장 + 통과대역 지표 자동 계산
사용법:
  1) 이 스크립트를 결과 hdf5 파일들과 같은 폴더에 둔다
  2) THROUGH / DROP 에 실제 flux 모니터 이름 입력
     (처음 실행 시 모니터 이름 목록이 출력되니 그걸 보고 고쳐도 됨)
  3) python export_crow_csv.py
링-링 gap과 격자 설정은 각 파일의 구조 정보에서 자동으로 읽는다.
"""
import glob
import os
import numpy as np
import pandas as pd
import tidy3d as td

THROUGH = "through"  # through 포트 flux 모니터 이름
DROP = "drop"        # drop 포트 flux 모니터 이름

# 통과대역 지표를 계산할 범위 (1550 nm 근처 공진만)
LAM_MIN, LAM_MAX = 1.546, 1.552

os.makedirs("data", exist_ok=True)


def ring_gap_from_sim(sim):
    """원기둥(링) 구조의 중심 y와 outer radius로 링-링 gap 계산"""
    cyl = [(s.geometry.center[1], s.geometry.radius) for s in sim.structures
           if isinstance(s.geometry, td.Cylinder)]
    if not cyl:
        return np.nan
    ys = np.array([c[0] for c in cyl])
    r_out = max(c[1] for c in cyl)
    return round((ys.max() - ys.min()) - 2 * r_out, 4)


def grid_from_sim(sim):
    try:
        return sim.grid_spec.auto_grid.min_steps_per_wvl
    except AttributeError:
        return np.nan


def passband_metrics(lam, d_db):
    """drop 스펙트럼(dB)에서 피크 2개, 최저점, 리플, 3 dB 대역폭 계산"""
    m = (lam >= LAM_MIN) & (lam <= LAM_MAX)
    lam, d = lam[m], d_db[m]

    # 국소 최대값 중 상위 2개를 피크로 사용
    idx = np.where((d[1:-1] > d[:-2]) & (d[1:-1] > d[2:]))[0] + 1
    top = idx[np.argsort(d[idx])[-2:]] if len(idx) >= 2 else idx
    top = np.sort(top)

    res = {}
    if len(top) == 2 and d[top[0]:top[1] + 1].min() < min(d[top]) - 0.01:
        seg = slice(top[0], top[1] + 1)
        i_dip = top[0] + np.argmin(d[seg])
        peak_avg = d[top].mean()
        res.update(
            peak1_um=lam[top[0]], peak1_dB=d[top[0]],
            peak2_um=lam[top[1]], peak2_dB=d[top[1]],
            split_nm=(lam[top[1]] - lam[top[0]]) * 1e3,
            dip_um=lam[i_dip], dip_dB=d[i_dip],
            ripple_dB=peak_avg - d[i_dip],
        )
        ref = d[top].max()
    else:  # 단일 피크
        i_pk = np.argmax(d)
        res.update(peak1_um=lam[i_pk], peak1_dB=d[i_pk], ripple_dB=0.0)
        ref = d[i_pk]

    # 3 dB 대역폭: 최대값 대비 -3 dB 교차점 (선형 보간)
    lvl = ref - 3.0
    above = np.where(d >= lvl)[0]
    i0, i1 = above[0], above[-1]
    lam_lo = np.interp(lvl, [d[i0 - 1], d[i0]], [lam[i0 - 1], lam[i0]])
    lam_hi = np.interp(lvl, [d[i1 + 1], d[i1]], [lam[i1 + 1], lam[i1]])
    res.update(bw3dB_lo_um=lam_lo, bw3dB_hi_um=lam_hi,
               bw3dB_nm=(lam_hi - lam_lo) * 1e3,
               center_um=(lam_lo + lam_hi) / 2)
    return res


summary = []
for h5 in sorted(glob.glob("*.hdf5")):
    try:
        sim_data = td.SimulationData.from_file(h5)
    except Exception as e:
        print(f"[skip] {h5}: 결과 파일이 아님 ({type(e).__name__})")
        continue

    sim = sim_data.simulation
    gap = ring_gap_from_sim(sim)
    grid = grid_from_sim(sim)
    names = [mon.name for mon in sim.monitors]
    print(f"[{h5}] ring_gap={gap} um, grid={grid}, monitors: {names}")

    out = {}
    for label, name in (("through", THROUGH), ("drop", DROP)):
        flux = sim_data[name].flux
        lam = td.C_0 / flux.coords["f"].values
        p = np.abs(flux.values)
        out["wavelength_um"] = lam
        out[f"{label}_lin"] = p
        out[f"{label}_dB"] = 10 * np.log10(p)

    df = pd.DataFrame(out).sort_values("wavelength_um").reset_index(drop=True)
    tag = f"gap{gap:.2f}_grid{grid}"
    df.to_csv(f"data/crow_{tag}.csv", index=False)

    met = passband_metrics(df["wavelength_um"].values, df["drop_dB"].values)
    met.update(file=h5, gap_um=gap, grid=grid)
    summary.append(met)
    print(f"   -> data/crow_{tag}.csv, ripple={met['ripple_dB']:.2f} dB, "
          f"BW3dB={met['bw3dB_nm']:.3f} nm")

cols = ["file", "gap_um", "grid", "peak1_um", "peak1_dB", "peak2_um",
        "peak2_dB", "split_nm", "dip_um", "dip_dB", "ripple_dB",
        "bw3dB_lo_um", "bw3dB_hi_um", "bw3dB_nm", "center_um"]
sm = pd.DataFrame(summary).reindex(columns=cols).sort_values(["grid", "gap_um"])
sm.to_csv("data/crow_summary.csv", index=False)
print(sm.to_string(index=False))
