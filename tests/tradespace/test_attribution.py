"""
Unit tests for the attribution of variance to uncertain parameters.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Choice,
    DesignSpace,
    Infeasible,
    Integer,
    ResultStore,
    Study,
    Uncertain,
    attribute_variance,
)

SPACE = DesignSpace(
    x=Integer(1, 2),
    u1=Uncertain.normal(0, 1),
    u2=Uncertain.normal(0, 1),
    s=Uncertain.seed(),
)


def _evaluate(d):
    """Metrics with known shares: y = x (u1 + 2 u2) (0.2 and 0.8 of u1 and u2,
    without interactions), z = u1 u2 (only their interaction), v (a draw of a
    seed), and w (constant)."""
    return {
        "y": d.x * (d.u1 + 2 * d.u2),
        "z": d.u1 * d.u2,
        "v": float(np.random.default_rng(d.s).normal()),
        "w": 1.0,
    }


def _share(attribution, metric, parameter, share="total", design=None):
    shares = attribution.shares
    rows = shares[(shares.metric == metric) & (shares.parameter == parameter)]
    if design is not None:
        rows = rows[rows.design == design]
    return rows[share].to_numpy()


class TestAttributeVariance(unittest.TestCase):
    """
    Unit tests for attributing variance to uncertain parameters.
    """

    @classmethod
    def setUpClass(cls):
        cls.study = Study(SPACE, lambda d: d, _evaluate, replicates=2)
        cls.designs = SPACE.enumerate()
        cls.attribution = attribute_variance(
            cls.study,
            cls.designs,
            ["y", "z", "v", "w"],
            replicates=300,
            n_resamples=200,
        )

    def test_known_shares(self):
        """The shares of additive effects and interactions are estimated (within
        three standard errors), and those of parameters without effects are 0."""
        a = self.attribution
        expected = {
            ("y", "u1"): (0.2, 0.2),
            ("y", "u2"): (0.8, 0.8),
            ("z", "u1"): (0.0, 1.0),
            ("z", "u2"): (0.0, 1.0),
            ("v", "s"): (1.0, 1.0),
        }
        for (metric, parameter), values in expected.items():
            for share, value in zip(("first", "total"), values):
                with self.subTest(metric=metric, parameter=parameter, share=share):
                    estimate = _share(a, metric, parameter, share)
                    error = _share(a, metric, parameter, f"{share} se")
                    self.assertTrue((np.abs(estimate - value) <= 3 * error).all())
        np.testing.assert_array_equal(_share(a, "y", "s"), 0.0)
        np.testing.assert_array_equal(_share(a, "z", "s"), 0.0)
        np.testing.assert_array_equal(_share(a, "v", "u1"), 0.0)
        np.testing.assert_array_equal(_share(a, "v", "s", "first"), 1.0)
        # with common random numbers, the shares of y do not depend on x
        np.testing.assert_allclose(*_share(a, "y", "u1").reshape(2, 1))

    def test_standard_errors(self):
        """The standard errors are those of the shares' estimates."""
        shares = self.attribution.shares
        y = shares[(shares.metric == "y") & (shares.parameter != "s")]
        self.assertTrue((y["total se"] > 0).all() and (y["total se"] < 0.1).all())
        # a share of 0 or 1 that every resample gives has no error
        v = shares[(shares.metric == "v") & (shares.parameter == "s")]
        np.testing.assert_array_equal(v["first se"], 0.0)

    def test_constant_metric(self):
        """The shares of a metric that does not vary are undefined."""
        shares = self.attribution.shares
        w = shares[shares.metric == "w"]
        self.assertTrue((w["variance"] == 0).all())
        self.assertTrue(
            w[["first", "total", "first se", "total se"]].isna().all().all()
        )

    def test_results(self):
        """The results have the plain and the mixed replicates."""
        results = self.attribution.results
        plain = results[results.redrawn.isna()]
        self.assertEqual(sorted(plain.replicate.unique()), list(range(600)))
        for name in ("u1", "u2", "s"):
            mixed = results[results.redrawn == name]
            self.assertEqual(len(mixed), 2 * 300)
            self.assertEqual(sorted(mixed.replicate.unique()), list(range(300)))
        # a mixed replicate draws its parameter as in B, and others as in A
        a = plain.set_index(["design", "replicate"])
        mixed = results[results.redrawn == "u1"].set_index(["design", "replicate"])
        shifted = [(d, r + 300) for d, r in mixed.index]
        np.testing.assert_array_equal(mixed.u1, a.loc[shifted].u1)
        np.testing.assert_array_equal(mixed.u2, a.loc[mixed.index].u2)

    def test_summary(self):
        """The summary pools the designs' shares, weighted by their variances."""
        a = self.attribution
        summary = a.summary()
        metrics = ["y", "z", "v", "w"]
        self.assertEqual(
            list(summary.columns), [*metrics, *(f"{m} se" for m in metrics)]
        )
        self.assertEqual(list(summary.index), ["u1", "u2", "s"])
        rows = a.shares[(a.shares.metric == "y") & (a.shares.parameter == "u2")]
        expected = (rows.total * rows.variance).sum() / rows.variance.sum()
        self.assertAlmostEqual(summary.loc["u2", "y"], expected)
        self.assertTrue(summary[["w", "w se"]].isna().all().all())
        # with common random numbers, the shares of y of both designs are the
        # same in every resample, and so is their pooled share
        np.testing.assert_allclose(summary.loc["u2", "y se"], rows["total se"])
        self.assertTrue((summary.loc[["u1", "u2"], ["y se", "z se"]] > 0).all().all())
        np.testing.assert_array_equal(summary.loc["s", ["y se", "z se"]], 0.0)
        first = a.summary("first")
        self.assertAlmostEqual(first.loc["s", "v"], 1.0)
        self.assertAlmostEqual(first.loc["s", "v se"], 0.0)
        self.assertEqual(a.samples["total"].shape, (len(a.shares), 200))
        with self.assertRaises(ValueError):
            a.summary("second")

    def test_groups(self):
        """Groups of parameters are redrawn together."""
        a = attribute_variance(
            self.study,
            self.designs[:1],
            ["y"],
            parameters={"normal": ["u1", "u2"], "seed": "s"},
            replicates=50,
            n_resamples=50,
        )
        self.assertEqual(list(a.shares.parameter), ["normal", "seed"])
        np.testing.assert_allclose(_share(a, "y", "normal", "first"), 1.0)
        np.testing.assert_allclose(_share(a, "y", "normal", "total"), 1.0, atol=0.3)
        np.testing.assert_array_equal(_share(a, "y", "seed", "total"), 0.0)
        named = attribute_variance(
            self.study, self.designs[:1], ["y"], parameters=["u2"], replicates=50
        )
        self.assertEqual(list(named.shares.parameter), ["u2"])
        self.assertEqual(set(named.results.redrawn.dropna()), {"u2"})

    def test_validation(self):
        """Unknown, repeated, or missing parameters and too few replicates raise."""
        for parameters in (["u3"], {"a": ["u1"], "b": ["u1", "u2"]}, {"a": []}, {}):
            with self.subTest(parameters=parameters):
                with self.assertRaises(ValueError):
                    attribute_variance(
                        self.study, self.designs, ["y"], parameters, replicates=2
                    )
        with self.assertRaises(ValueError):
            attribute_variance(self.study, self.designs, ["y"], replicates=1)
        with self.assertRaises(ValueError):
            attribute_variance(self.study, self.designs, ["q"], replicates=2)
        certain = Study(DesignSpace(x=Integer(1, 2)), lambda d: d, lambda d: {"y": d.x})
        with self.assertRaises(ValueError):
            attribute_variance(certain, certain.space.enumerate(), ["y"])

    def test_inactive_branches(self):
        """Designs that do not draw a parameter are not evaluated for it."""
        space = DesignSpace(
            layout=Choice(
                a=DesignSpace(phase=Uncertain.uniform(0, 1)),
                b=DesignSpace(size=Integer(1, 2)),
            ),
            noise=Uncertain.uniform(0, 1),
        )

        def evaluate(d):
            phase = d.layout.phase if d.layout.choice == "a" else 0.0
            return {"y": phase + d.noise}

        study = Study(space, lambda d: d, evaluate)
        designs = space.enumerate()
        a = attribute_variance(study, designs, ["y"], replicates=20, n_resamples=20)
        redrawn = a.results[a.results.redrawn == "layout.a.phase"]
        self.assertEqual(set(redrawn.design), {designs[0].key})
        shares = a.shares.set_index(["design", "parameter"])
        for design in designs[1:]:
            self.assertTrue(
                np.isnan(shares.loc[(design.key, "layout.a.phase"), "total"])
            )
            self.assertAlmostEqual(
                shares.loc[(design.key, "noise"), "total"], 1.0, delta=0.5
            )

    def test_infeasible(self):
        """Infeasible evaluations are left out."""

        def build(d):
            if d.u1 > 1:
                raise Infeasible("too large")
            return d

        study = Study(SPACE, build, _evaluate)
        a = attribute_variance(
            study, self.designs[:1], ["y"], parameters=["u1", "u2"], replicates=100
        )
        self.assertFalse(a.results.feasible.all())
        self.assertTrue(a.shares[["first", "total"]].notna().all().all())

        def reject(d):
            raise Infeasible("never")

        never = Study(SPACE, reject, _evaluate)
        with self.assertWarns(UserWarning):
            none = attribute_variance(never, self.designs, ["y"], replicates=2)
        self.assertEqual(len(none.shares), 0)
        self.assertEqual(len(none.summary()), 0)
        with self.assertWarns(UserWarning):
            separated = attribute_variance(
                never, self.designs, ["y"], replicates=2, by=["sensor"]
            )
        self.assertEqual(len(separated.shares), 0)

    def test_by(self):
        """The records of an evaluation are separated by `by`."""

        def evaluate(d):
            return pd.DataFrame(
                {"y": [d.u1, d.u2]}, index=pd.Index(["one", "two"], name="sensor")
            )

        study = Study(SPACE, lambda d: d, evaluate)
        a = attribute_variance(
            study,
            self.designs[:1],
            ["y"],
            parameters=["u1", "u2"],
            replicates=20,
            by=["sensor"],
        )
        shares = a.shares.set_index(["sensor", "parameter"])["total"]
        self.assertAlmostEqual(shares["one", "u1"], 1.0, delta=0.5)
        self.assertEqual(shares["one", "u2"], 0.0)
        self.assertEqual(shares["two", "u1"], 0.0)
        self.assertEqual(list(a.summary().index.names), ["sensor", "parameter"])
        with self.assertRaises(ValueError):
            attribute_variance(study, self.designs[:1], ["y"], replicates=2)

    def test_store(self):
        """The evaluations are stored, and taken from the store again."""
        self.calls = []

        def evaluate(d):
            self.calls.append(d.replicate)
            return _evaluate(d)

        with tempfile.TemporaryDirectory() as directory:
            store = ResultStore(Path(directory) / "results.sqlite", "attribution")
            study = Study(SPACE, lambda d: d, evaluate, replicates=3, store=store)
            first = attribute_variance(study, self.designs, ["y"], replicates=3)
            self.assertEqual(len(self.calls), 2 * (2 * 3 + 3 * 3))
            self.calls.clear()
            again = attribute_variance(study, self.designs, ["y"], replicates=3)
            self.assertEqual(self.calls, [])
            pd.testing.assert_frame_equal(first.shares, again.shares)


def _evaluate_pairs(d):
    """Metrics of differences with known shares: y = x (u1 + 2 u2), whose
    difference between the designs is that of one design, and c = x u1 + u2,
    whose difference is u1 alone (u2 cancels)."""
    return {"y": d.x * (d.u1 + 2 * d.u2), "c": d.x * d.u1 + d.u2}


class TestDifferences(unittest.TestCase):
    """
    Unit tests for attributing the variance of differences between designs.
    """

    @classmethod
    def setUpClass(cls):
        cls.designs = SPACE.enumerate()
        cls.study = Study(SPACE, lambda d: d, _evaluate_pairs)
        cls.attribution = attribute_variance(
            cls.study,
            cls.designs,
            ["y", "c"],
            parameters=["u1", "u2"],
            replicates=200,
            n_resamples=200,
        )
        cls.differences = cls.attribution.differences()

    def test_cancelled(self):
        """Parameters that act alike on both designs cancel in differences."""
        d = self.differences
        self.assertEqual(list(d.shares.columns[:2]), ["design", "other"])
        self.assertEqual(
            (d.shares.design.iloc[0], d.shares.other.iloc[0]),
            (self.designs[0].key, self.designs[1].key),
        )
        # c differs by u1 alone: u2 explains none of its variance in any resample
        np.testing.assert_allclose(_share(d, "c", "u2"), 0.0, atol=1e-12)
        np.testing.assert_allclose(_share(d, "c", "u1", "first"), 1.0)
        estimate, error = _share(d, "c", "u1"), _share(d, "c", "u1", "total se")
        self.assertTrue((np.abs(estimate - 1) <= 3 * error).all())
        self.assertAlmostEqual(d.summary().loc["u2", "c se"], 0.0)
        # u2 spreads c of both designs, but not their difference
        self.assertTrue((_share(self.attribution, "c", "u2") > 0.1).all())

    def test_differences_of_one_design(self):
        """The difference of y is that of a design, whose shares it has."""
        d, a = self.differences, self.attribution
        for share in ("first", "total", "first se", "total se"):
            for parameter in ("u1", "u2"):
                np.testing.assert_allclose(
                    _share(d, "y", parameter, share),
                    _share(a, "y", parameter, share, self.designs[0].key),
                )

    def test_results(self):
        """The results are the differences in each replicate."""
        results, differences = self.attribution.results, self.differences.results
        plain = results[results.redrawn.isna()].set_index(["design", "replicate"])
        first, other = self.designs[0].key, self.designs[1].key
        expected = plain.loc[first, "c"] - plain.loc[other, "c"]
        found = differences[differences.redrawn.isna()].set_index("replicate")["c"]
        np.testing.assert_allclose(found.sort_index(), expected.sort_index())
        self.assertEqual(len(differences), 2 * 200 + 2 * 200)

    def test_pairs(self):
        """Pairs are given by designs or keys, and reversed pairs have the
        same shares."""
        first, other = self.designs
        reversed_ = self.attribution.differences([(other, first.key)])
        columns = ["variance", "first", "total", "first se", "total se"]
        np.testing.assert_allclose(
            reversed_.shares[columns], self.differences.shares[columns]
        )
        self.assertEqual(reversed_.shares.design.iloc[0], other.key)
        for pairs in ([(first, first)], [(first.key, "unknown")]):
            with self.subTest(pairs=pairs):
                with self.assertRaises(ValueError):
                    self.attribution.differences(pairs)
        with self.assertRaises(ValueError):
            self.differences.differences()

    def test_inactive_branches(self):
        """A design that does not draw a parameter is as in A."""
        space = DesignSpace(
            layout=Choice(
                a=DesignSpace(phase=Uncertain.uniform(0, 1)),
                b=DesignSpace(size=Integer(1, 2)),
            ),
            noise=Uncertain.uniform(0, 1),
        )

        def evaluate(d):
            if d.layout.choice == "a":
                return {"y": d.layout.phase + d.noise}
            return {"y": d.layout.size * d.noise}

        study = Study(space, lambda d: d, evaluate)
        designs = space.enumerate()
        a = attribute_variance(study, designs, ["y"], replicates=20, n_resamples=20)
        d = a.differences().shares.set_index(["design", "other", "parameter"])
        branch_a, b1, b2 = (design.key for design in designs)
        # a less b1 is the phase alone
        self.assertAlmostEqual(d.loc[(branch_a, b1, "noise"), "total"], 0.0)
        self.assertAlmostEqual(d.loc[(branch_a, b1, "layout.a.phase"), "first"], 1.0)
        # neither b1 nor b2 draws the phase
        self.assertTrue(np.isnan(d.loc[(b1, b2, "layout.a.phase"), "total"]))
        self.assertAlmostEqual(d.loc[(b1, b2, "noise"), "first"], 1.0)

    def test_independent_draws(self):
        """Without common random numbers, no parameter cancels."""
        study = Study(SPACE, lambda d: d, _evaluate_pairs, common_random_numbers=False)
        a = attribute_variance(
            study, self.designs, ["c"], parameters=["u1", "u2"], replicates=200
        )
        d = a.differences()
        # c of the designs differs by u1 of both and u2 of both: 5 and 2 of 7
        for parameter, value in (("u1", 5 / 7), ("u2", 2 / 7)):
            estimate, error = _share(d, "c", parameter), _share(
                d, "c", parameter, "total se"
            )
            self.assertTrue((np.abs(estimate - value) <= 3 * error).all())


if __name__ == "__main__":
    unittest.main()
