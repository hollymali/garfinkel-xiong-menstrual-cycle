"""Numerical smoke tests for the Meeting 10 sleep/SCN edition."""

import numpy as np

import MODEL


def main():
    boundary_zt = np.array([0.0, 1.0, 13.0, 15.999, 16.0, 23.999])
    expected_sleep = np.array([1.0, 1.0, 1.0, 1.0, 1.57, 1.57])
    np.testing.assert_allclose(
        MODEL.sleep_period_multiplier(boundary_zt),
        expected_sleep,
    )

    scn_at_peak_and_trough = MODEL.scn_drive(np.array([1.0, 13.0]))
    np.testing.assert_allclose(scn_at_peak_and_trough, np.array([2.0, 1.0]))

    trajectory = MODEL.run(output_h=24.0)
    scn = trajectory["SCN_drive"]
    assert scn.min() >= 1.0 - 1e-12
    assert scn.max() <= 2.0 + 1e-12

    # Count counter-crossings in one retained day. This checks that the sleep
    # multiplier is affecting timing rather than only being returned as data.
    phi = trajectory["phi"]
    pulse_indices = np.flatnonzero(np.diff(np.floor(phi)) > 0.0) + 1
    pulse_zt = trajectory["zt_h"][pulse_indices]
    sleep_pulses = np.count_nonzero((pulse_zt >= 16.0) & (pulse_zt < 24.0))
    wake_pulses = len(pulse_zt) - sleep_pulses
    sleep_mask = (trajectory["zt_h"] >= 16.0) & (trajectory["zt_h"] < 24.0)
    wake_mask = ~sleep_mask

    print("clock boundary checks: passed")
    print(f"SCN min/max: {scn.min():.3f}/{scn.max():.3f}")
    print(f"one-day pulse count wake/sleep: {wake_pulses}/{sleep_pulses}")
    print(
        "mean E2 wake/sleep: "
        f"{trajectory['E2_pg_ml'][wake_mask].mean():.3f}/"
        f"{trajectory['E2_pg_ml'][sleep_mask].mean():.3f} pg/mL"
    )
    print(
        "mean LH wake/sleep: "
        f"{trajectory['LH_IU_L'][wake_mask].mean():.3f}/"
        f"{trajectory['LH_IU_L'][sleep_mask].mean():.3f} IU/L"
    )
    print(
        "E2 min/mean/max: "
        f"{trajectory['E2_pg_ml'].min():.3f}/"
        f"{trajectory['E2_pg_ml'].mean():.3f}/"
        f"{trajectory['E2_pg_ml'].max():.3f} pg/mL"
    )


if __name__ == "__main__":
    main()
