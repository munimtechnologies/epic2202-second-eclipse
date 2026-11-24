# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Stdlib tests for the EPIC 220208795 second-eclipse test (hand-traced values)."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ecl  # noqa: E402
import io_ztf  # noqa: E402


class TestTemplate(unittest.TestCase):
    def test_values(self):
        self.assertAlmostEqual(ecl.template_flux(-0.0250), 0.7166, places=6)
        self.assertAlmostEqual(ecl.template_flux(-0.0352), 0.5 * (0.7447 + 0.7166), places=6)
        self.assertEqual(ecl.template_flux(-1.0), 1.0)
        self.assertEqual(ecl.template_flux(0.5), 1.0)
        self.assertAlmostEqual(ecl.template_depth(-0.0250), 0.2834, places=6)

    def test_core_edges(self):
        # depth just below 0.15 at both edges, every template point with depth >= 0.15 inside
        self.assertLess(ecl.template_depth(ecl.CORE_LO), ecl.CORE_DEPTH)
        self.assertLess(ecl.template_depth(ecl.CORE_HI), ecl.CORE_DEPTH)
        self.assertAlmostEqual(ecl.template_flux(ecl.CORE_LO), 0.852567, places=5)
        for u, f in ecl.TEMPLATE:
            if 1 - f >= ecl.CORE_DEPTH:
                self.assertTrue(ecl.CORE_LO < u < ecl.CORE_HI)
        self.assertGreaterEqual(ecl.template_depth(-0.0659), ecl.CORE_DEPTH)
        self.assertGreaterEqual(ecl.template_depth(0.0363), ecl.CORE_DEPTH)

    def test_dates(self):
        self.assertEqual(ecl.bjd_to_date(2451545.0), "2000-01-01 12:00")
        self.assertEqual(ecl.bjd_to_date(2460942.764), "2025-09-24 06:20")


def _pair(T, r=1.0, s=0.01, src="ztf"):
    return [(T - 0.025, r, s, src), (T - 0.0046, r, s, src)]


class TestVeto(unittest.TestCase):
    def test_two_normal_points_veto(self):
        T = ecl.T_REF + 300.0
        pts = ecl.Points(_pair(T))
        v = ecl.vetoes(pts, T)
        self.assertEqual(len(v), 2)

    def test_dipped_points_do_not_veto(self):
        T = ecl.T_REF + 300.0
        pts = ecl.Points(_pair(T, r=0.75))
        self.assertEqual(len(pts.vt), 0)
        self.assertEqual(ecl.vetoes(pts, T), [])

    def test_noisy_point_does_not_veto(self):
        # r - f = 1.0 - 0.7166 = 0.2834 < 4 sigma when sigma = 0.08
        T = ecl.T_REF + 300.0
        pts = ecl.Points([(T - 0.025, 1.0, 0.08, "ztf")])
        self.assertEqual(ecl.vetoes(pts, T), [])

    def test_outside_core_does_not_veto(self):
        T = ecl.T_REF + 300.0
        pts = ecl.Points([(T + 0.15, 1.0, 0.01, "ztf"), (T - 0.15, 1.0, 0.01, "ztf")])
        self.assertEqual(ecl.vetoes(pts, T), [])

    def test_source_thresholds(self):
        pts = ecl.Points([(1.0, 0.96, 0.01, "tess"), (2.0, 0.96, 0.01, "ztf")])
        self.assertEqual(pts.vsrc, ["ztf"])

    def test_trial_and_scan(self):
        T = ecl.T_REF + 300.0
        pts = ecl.Points(_pair(T))
        ex, inf = ecl.trial(300.0, pts, ecl.T_REF - 10, ecl.T_REF + 1000)
        self.assertTrue(ex)
        self.assertEqual(inf["n"], 1)
        self.assertEqual(inf["by"], ["ztf"])
        ex2, _ = ecl.trial(310.0, pts, ecl.T_REF - 10, ecl.T_REF + 1000)
        self.assertFalse(ex2)
        # one point only: excluded with min_veto 1, not with 2
        p1 = ecl.Points([(T - 0.025, 1.0, 0.01, "ztf")])
        self.assertFalse(ecl.trial(300.0, p1, ecl.T_REF - 10, ecl.T_REF + 1000)[0])
        self.assertTrue(ecl.trial(300.0, p1, ecl.T_REF - 10, ecl.T_REF + 1000, min_veto=1)[0])
        # P = 150 also predicts an eclipse at T (n = 2)
        self.assertTrue(ecl.trial(150.0, pts, ecl.T_REF - 10, ecl.T_REF + 1000)[0])

    def test_eclipse_refs(self):
        refs = ecl.eclipse_refs(100.0, ecl.T_REF - 250, ecl.T_REF + 250)
        self.assertEqual([n for n, _ in refs], [-2, -1, 1, 2])
        self.assertAlmostEqual(refs[0][1], ecl.T_REF - 200.0)

    def test_scan_runs(self):
        T = ecl.T_REF + 300.0
        pts = ecl.Points(_pair(T))
        grid = [299.8, 299.9, 300.0, 300.1, 300.2]
        # both points veto only for T in about [Tc - 0.057, Tc + 0.050] -> only P = 300.0 excluded
        n_ex, runs = ecl.scan(pts, grid, ecl.T_REF - 10, ecl.T_REF + 400)
        self.assertEqual(n_ex, 1)
        self.assertEqual(runs, [[299.8, 299.9, 2], [300.1, 300.2, 2]])
        fr = ecl.binned_fraction(grid, runs, [299.0, 300.0, 301.0])
        self.assertEqual(fr[0][2:4], (2, 2))
        self.assertEqual(fr[1][2:4], (3, 2))

    def test_period_grid(self):
        g = ecl.period_grid(60.0, 61.0, ecl.T_REF - 100, ecl.T_REF + 1000)
        self.assertEqual(g[0], 60.0)
        self.assertLessEqual(g[-1], 61.0)
        fac = 1 + ecl.GRID_TOL / 1000.0
        self.assertAlmostEqual(g[1] / g[0], fac, places=12)
        expected = int(math.floor(math.log(61.0 / 60.0) / math.log(fac))) + 1
        self.assertIn(len(g), (expected - 1, expected, expected + 1))


class TestEvents(unittest.TestCase):
    def test_dips_and_groups(self):
        pts = ecl.Points([(10.0, 0.8, 0.01, "ztf"), (10.1, 0.78, 0.01, "asassn"),
                          (50.0, 0.85, 0.01, "ztf"), (60.0, 0.95, 0.01, "ztf"),
                          (70.0, 0.85, 0.05, "ztf"), (80.0, 0.7, 0.01, "tess")])
        d = ecl.dip_points(pts)
        self.assertEqual([x[0] for x in d], [10.0, 10.1, 50.0])
        ev = ecl.group_events(d)
        self.assertEqual(len(ev), 2)
        self.assertTrue(ev[0]["strong"])
        self.assertEqual(ev[0]["n"], 2)
        self.assertEqual(ev[0]["sources"], ["asassn", "ztf"])
        self.assertAlmostEqual(ev[0]["depth_max"], 0.22)
        self.assertFalse(ev[1]["strong"])

    def test_periods_from_event(self):
        t_min = ecl.T_REF + ecl.T_MIN_OFFSET + 600.0
        ps = ecl.periods_from_event(t_min)
        self.assertEqual(ps[0][0], 1)
        self.assertAlmostEqual(ps[0][1], 600.0, places=6)
        self.assertIn(len(ps), (9, 10))
        self.assertTrue(any(n == 2 and abs(P - 300.0) < 1e-6 for n, P in ps))
        neg = ecl.periods_from_event(ecl.T_REF + ecl.T_MIN_OFFSET - 200.0)
        self.assertEqual(neg[0][0], -1)
        self.assertAlmostEqual(neg[0][1], 200.0, places=6)

    def test_verdict(self):
        self.assertEqual(ecl.candidate_verdict(1.0, True), "SUPPORTED")
        self.assertEqual(ecl.candidate_verdict(0.995, False), "EXCLUDED")
        self.assertEqual(ecl.candidate_verdict(0.98, False), "OPEN")


class TestSeries(unittest.TestCase):
    def _series(self):
        t = [k * 0.01 for k in range(1000)]
        r = [1.0 + 0.001 * (((k * 7) % 5) - 2) / 2.0 for k in range(1000)]
        return t, r

    def test_box_search_finds_injected(self):
        t, r = self._series()
        ri = [x * ecl.template_flux(ti - 5.0) for ti, x in zip(t, r)]
        best, dets = ecl.box_search(t, ri, 0.01)
        self.assertEqual(len(dets), 1)
        self.assertLess(abs(dets[0]["t_mid"] - 4.975), 0.08)
        self.assertGreater(dets[0]["depth"], 0.15)
        self.assertLess(abs(best["t_mid"] - 4.975), 0.08)

    def test_box_search_flat(self):
        t, r = self._series()
        best, dets = ecl.box_search(t, r, 0.01)
        self.assertEqual(dets, [])
        self.assertLess(best["depth"], 0.002)

    def test_normalise_and_bin(self):
        t, _ = self._series()
        f = [2.0] * len(t)
        r = ecl.normalise_binned(t, f)
        self.assertTrue(all(abs(x - 1.0) < 1e-12 for x in r))
        # 2026-11-24 fix: the first version used times 0.5 and 0.51 with width 0.1; in floating
        # point 0.5 // 0.1 == 4.0 but 0.51 // 0.1 == 5.0, so those two fell in different bins
        # (the code was right, the test input sat on a bin edge). Times now sit mid-bin.
        b = ecl.bin_series([0.0, 0.01, 0.02, 0.53, 0.55], [1.0, 2.0, 3.0, 4.0, 6.0], 0.1)
        self.assertEqual(len(b), 2)
        self.assertAlmostEqual(b[0][1], 2.0)
        self.assertEqual(b[0][2], 3)
        self.assertAlmostEqual(b[1][1], 5.0)
        self.assertEqual(b[1][2], 2)

    def test_robust_scatter(self):
        self.assertAlmostEqual(ecl.robust_scatter([1.0, 2.0, 3.0, 4.0, 5.0]), 1.4826, places=6)


class TestZtfParse(unittest.TestCase):
    def test_parse(self):
        head = "oid,expid,hjd,mjd,mag,magerr,catflags,filtercode\n"
        rows = ["A,%d,%.5f,%.5f,15.0,0.01,0,zg" % (k, 2458300.0 + k, 58299.5 + k) for k in range(11)]
        rows.append("A,99,2458400.00000,58399.50000,15.3,0.01,0,zg")
        rows.append("A,98,2458401.00000,58400.50000,15.0,0.01,4,zg")
        rows.append("B,97,2458402.00000,58401.50000,14.0,0.01,0,zr")
        pts, meta = io_ztf.parse(head + "\n".join(rows) + "\n")
        self.assertEqual(len(pts), 12)
        self.assertEqual(meta["n_raw"], 14)
        self.assertFalse(meta["oids"]["B"]["used"])
        rmin = min(p[1] for p in pts)
        self.assertAlmostEqual(rmin, 10 ** (-0.4 * 0.3), places=9)
        self.assertAlmostEqual(pts[0][2], 0.00921, places=9)
        self.assertEqual(pts[0][3], "ztf")


class TestAsassnNormalise(unittest.TestCase):
    def test_normalise(self):
        import asassn_io
        by = {"g": [(2459000.0 + k, 10.0, 0.1, "bA") for k in range(25)]
              + [(2459100.0 + k, 5.0, 0.05, "bB") for k in range(5)], "V": []}
        pts, meta = asassn_io.normalise(by)
        self.assertEqual(len(pts), 25)
        self.assertTrue(meta["g:bA"]["used"])
        self.assertFalse(meta["g:bB"]["used"])
        self.assertAlmostEqual(pts[0][1], 1.0)
        self.assertAlmostEqual(pts[0][2], 0.01)


if __name__ == "__main__":
    unittest.main()
