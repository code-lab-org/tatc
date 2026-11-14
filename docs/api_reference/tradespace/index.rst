===========
Trade Space
===========

.. automodule:: tatc.tradespace
   :no-members:

The :doc:`../../examples/tradespace/index` shows each step on a study of the
observations of constellations, one step at a time.

Steps of a Trade Study
======================

A trade study takes a few steps, each of which takes the outcome of the one
before::

    from tatc.tradespace import (
        Criteria, DesignSpace, Integer, ResultStore, Study, Uncertain,
        aggregate, confirm, pareto_front, refine,
    )

    space = DesignSpace(planes=Integer(1, 4), per_plane=Integer(1, 4), seed=Uncertain.seed())
    study = Study(space, build, evaluate, replicates=5, store=ResultStore("results.sqlite"))
    results = study.run()  # a row per design and replicate

    criteria = Criteria({"satellites": "min", "coverage": "max"})
    summary = aggregate(results, criteria)  # a row per design
    front = summary[pareto_front(summary, criteria)]

    refinement = refine(study, criteria, results)  # replicates where efficiency is uncertain
    confirmation = confirm(study, criteria, refinement.results)  # the front, with new replicates

1. **Define** a design space (:doc:`define`): design variables, choices and
   repetitions, constraints, and uncertain parameters, which each replicate
   of a design draws.
2. **Evaluate** designs with a study (:doc:`evaluate`): a function builds
   each realized design and another evaluates it, in parallel, with a store
   of the evaluations. An `ObservationStudy` evaluates the observations of
   constellations, and their deployment, costs, and failures
   (:doc:`missions`).
3. **Judge** designs by criteria (:doc:`judge`): objectives and
   requirements on statistics of their replicates, summarized with a row per
   design, and their Pareto efficiency.
4. **Settle** the uncertainty of the replicates, or **search** spaces too
   large to enumerate (:doc:`search`). Each of these steps takes a study and
   criteria, and returns its `results`, the `summary` of its designs, and the
   `designs` themselves, for the next step.

Results and Summaries
=====================

A study gives **results**: a row per evaluation (a design in a replicate,
and a scenario), with the key of the design (`design`), the `replicate`, the
values of the design by path (see `Design.flat`), whether it is `feasible`,
and its metrics. A **summary** has a row per design (and per group of the
criteria's `by`): statistics over its replicates, the number of
`replicates`, and the values of the design. Functions take one or the
other:

.. list-table::
   :header-rows: 1

   * - Take results (a row per replicate)
     - Take summaries (a row per design)
   * - `aggregate`, `summarize`, `pareto_probability`,
       `hypervolume_progress`, `robustness`, `refine`, `confirm`,
       `attribute_variance`, `plot_means`
     - `pareto_front`, `pareto_ranks`, `hypervolume`, `sensitivity`,
       `main_effects`, `plot_designs`

The functions of summaries raise an error if a design has several rows, as
in the results of a study: summarize them first.

Columns
-------

Columns named after a metric or a statistic take a suffix after a space,
and columns that TAT-C defines are named in snake case:

* statistics over the replicates of a design: named by the criteria (such
  as `"cost p70"` for `Quantile("cost", 0.7)`), the means of metrics by
  their names, and `<metric> std` and `<metric> se` (by `summarize`, and the
  standard error of each objective by `pareto_probability`).
* statistics within a replicate, over the mission of a design (see
  `ObservationStudy`): `<metric> mission` (the mean over the phases of the
  mission), `<metric> time` (when the phases first meet a requirement), and,
  over the histories of its failures, `<column> histories p<percent>`,
  `<column> histories std`, and `P(<requirement>)`.
* columns of judgments: `violation`, `meets_requirements`, `pareto`,
  `pareto_rank`, `pareto_probability`, `requirements_probability`,
  `pareto_trace`, `max_regret`, and the `<scenario> regret` of each scenario.

Designs and Their Variants
==========================

Results keep the values of designs, from which `DesignSpace.designs_of`
rebuilds them, so that the steps of a study take designs as `Design`
objects, their keys, or rows of a summary or of results::

    confirmation = confirm(study, criteria, results, summary[summary.pareto])

A variant of a design space adds or replaces children
(`DesignSpace.extend`), and its designs take the decisions of a design of
the other (`DesignSpace.design`). A variant of a study replaces some of its
arguments (`Study.replace`)::

    policies = space.extend(policy=policy_space(orbit_spares=Integer(0, 1)))
    designs = [policies.design(d, policy=p) for d in front_designs for p in policy_designs]
    results = study.replace(space=policies, replicates=20).run(designs)

The parts of an `ObservationStudy` (such as its targets and mission plan)
refer to values of designs by a `Reference` (see
:py:data:`tatc.tradespace.observations.Reference`): a name or path of a
value, or a function of the realized design.

Stores of Results
=================

A `ResultStore` keeps the evaluations of studies, keyed by the fingerprint
of each study's definition (see `Study.fingerprint`), so that a study that is
run again takes its evaluations from the store, and variants of a study
(other settings, design spaces, or scenarios) share a store without mixing
their evaluations. A namespace of the store identifies only what
fingerprints cannot see (such as data that the functions read from files). A
study warns if the store has evaluations of the same study with other code
(of its functions, of TAT-C, or of Python), which it does not reuse
(`ResultStore.prune` removes them).

.. toctree::
  :maxdepth: 1

  define
  evaluate
  judge
  search
  missions
