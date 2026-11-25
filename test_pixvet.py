# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
import math
import unittest

import pixvet


class TestStats(unittest.TestCase):
    def test_median(self):
        self.assertEqual(pixvet.median([3, 1, 2]), 2.0)
        self.assertEqual(pixvet.median([4, 1, 3, 2]), 2.5)
        self.assertTrue(math.isnan(pixvet.median([])))
        self.assertEqual(pixvet.median([1.0, float("nan"), 3.0]), 2.0)

    def test_percentile(self):
        self.assertAlmostEqual(pixvet.percentile([1, 2, 3, 4, 5], 60), 3.4, places=9)
        self.assertAlmostEqual(pixvet.percentile([5, 1], 50), 3.0, places=9)

    def test_pearson(self):
        self.assertAlmostEqual(pixvet.pearson([1, 2, 3], [2, 4, 6]), 1.0, places=9)
        self.assertAlmostEqual(pixvet.pearson([1, 2, 3], [3, 2, 1]), -1.0, places=9)
        self.assertTrue(math.isnan(pixvet.pearson([1, 1, 1], [1, 2, 3])))

    def test_depth(self):
        self.assertAlmostEqual(pixvet.depth([0.9, 0.9, 0.9], [1.0, 1.0, 1.0]), 0.1, places=9)


class TestWindows(unittest.TestCase):
    def test_window_idx(self):
        ins, outs = pixvet.window_idx([0.0, 0.05, 0.1, 0.2, 0.4, 0.6], 0.0, 0.12)
        self.assertEqual(ins, [0, 1])
        self.assertEqual(outs, [3, 4])


def grid(n, val=0.0):
    return [[val for _ in range(n)] for _ in range(n)]


class TestImages(unittest.TestCase):
    def test_find_target(self):
        img = grid(5)
        img[2][3] = 10.0
        img[0][0] = 100.0
        self.assertEqual(pixvet.find_target(img), (2, 3))

    def test_bkg_mask(self):
        img = [[float(y * 11 + x) for x in range(11)] for y in range(11)]
        m = pixvet.bkg_mask(img, 5, 5, 5, 5)
        for (y, x) in m:
            self.assertLess(img[y][x], 72.0)
            self.assertGreaterEqual(math.hypot(y - 5, x - 5), 4.0)
        self.assertIn((0, 0), m)
        self.assertIn((5, 0), m)
        self.assertIn((6, 0), m)
        self.assertNotIn((5, 2), m)
        self.assertNotIn((6, 6), m)
        self.assertNotIn((9, 9), m)
        clipped = pixvet.bkg_mask(img, 0, 0, 0, 0)
        self.assertTrue(all(y <= 5 and x <= 5 for (y, x) in clipped))

    def test_ap_sum_and_box(self):
        img = grid(5, 1.0)
        self.assertEqual(pixvet.ap_sum(img, 2, 2), 9.0)
        self.assertEqual(pixvet.ap_sum(img, 0, 0), 4.0)
        b = pixvet.box(grid(3, 1.0), 0, 0, r=1)
        self.assertEqual(len(b), 9)
        self.assertEqual(sum(1 for v in b if math.isnan(v)), 5)

    def test_local_maxima(self):
        img = grid(11)
        for y in range(4, 7):
            for x in range(4, 7):
                img[y][x] = 10.0
        img[5][5] = 100.0
        img[2][8] = 50.0   # control, distance 4.24
        img[8][8] = 30.0   # control, fainter
        img[8][2] = 5.0    # too faint (3x3 sum 5 < 0.05 x 180)
        img[5][9] = 80.0   # too close to the edge
        self.assertEqual(pixvet.local_maxima(img, 5, 5), [(2, 8), (8, 8)])


class TestLabels(unittest.TestCase):
    def test_label_event(self):
        L = pixvet.label_event
        self.assertEqual(L(0.01, 0.9, 1.0, [0.0, 0.0]), "REPRO_FAIL")
        self.assertEqual(L(float("nan"), 0.9, 1.0, [0.0, 0.0]), "REPRO_FAIL")
        self.assertEqual(L(0.1, 0.9, 1.0, [0.0, 0.01, -0.01]), "ON_TARGET")
        self.assertEqual(L(0.1, 0.9, 1.0, []), "ON_TARGET")
        self.assertEqual(L(0.1, 0.9, 1.0, [0.08, 0.06]), "SYSTEMATIC")
        self.assertEqual(L(0.1, 0.1, 1.0, [0.0, 0.0]), "SYSTEMATIC")
        self.assertEqual(L(0.1, 0.9, 0.5, [0.0, 0.0]), "UNCLEAR")
        self.assertEqual(L(0.1, 0.5, 1.0, [0.0, 0.0]), "UNCLEAR")
        self.assertEqual(L(0.1, 0.9, 1.0, [0.04, 0.04]), "UNCLEAR")

    def test_gate_and_overall(self):
        self.assertEqual(pixvet.control_gate("ON_TARGET", "REPRO_FAIL"), "CONTROLS_OK")
        self.assertEqual(pixvet.control_gate("UNCLEAR", "REPRO_FAIL"), "CONTROL_FAIL")
        self.assertEqual(pixvet.control_gate("ON_TARGET", "SYSTEMATIC"), "CONTROL_FAIL")
        self.assertEqual(pixvet.overall("CONTROL_FAIL", {"a": "ON_TARGET"}), "UNINTERPRETED")
        self.assertEqual(pixvet.overall("CONTROLS_OK", {"a": "ON_TARGET", "b": "SYSTEMATIC"}),
                         "CANDIDATE_ECLIPSE")
        self.assertEqual(pixvet.overall("CONTROLS_OK", {"a": "SYSTEMATIC", "b": "SYSTEMATIC"}),
                         "ALL_SYSTEMATIC")
        self.assertEqual(pixvet.overall("CONTROLS_OK", {"a": "SYSTEMATIC", "b": "UNCLEAR"}), "MIXED")


if __name__ == "__main__":
    unittest.main()
