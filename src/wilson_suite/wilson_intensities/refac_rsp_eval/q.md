# Design memo — response-function evaluation boundaries

Working notes on where `wilson_intensities` ends and `wilson_derive` / `wilson_main` begin,
and what shape the evaluation stage should take. Written against the draft in `eval.py`.

---

## 1. The question

`eval.py` contains classes annotated `# symbolic` (`PropsCollection`, `FreqTermsCollection`,
`ResonanceMotif`) that hold derive-owned objects (`PolProp`, `VibDiffTerm`, `ResonanceCondition`)
but exist to serve evaluation. They look like they could belong to either package. Where do they go,
and how is that decided by a rule rather than by taste?

---

## 2. Diagnosis: the split is on the wrong axis

The file is partitioned `# symbolic` / `# numerical` — a classification by **what an object contains**.
The boundary-determining question is **what an object's value depends on**, because that determines
what goes stale, what can be cached, and what can be tested without a molecule.

There are three independent inputs:

| input | changes when |
|---|---|
| derived terms | order / experiment type changes |
| axis & experiment choice | the spectroscopy changes |
| molecular data | molecule / level of theory changes |

Two labels cannot partition three inputs. Hence the unresolvable feeling. The fix is a third stage,
already half-written as `parse_vibpert_term` → `RspEvalTerm`.

---

## 3. Ownership rule: identity decides

"Who calls it" is a weak test — it describes today's call sites, and a future caller can always
be imagined. The strong test:

> **A type belongs to the layer that defines its identity.** If `__eq__` / `__hash__` / the canonical
> form is chosen to serve a downstream storage or dispatch strategy, the type is downstream's,
> whatever it holds.

This is checkable from the class definition alone. And the codebase already settles the case.

**Derive's notion of identity**, for these same objects:

- `VibPerturbedTerm.sort(nm_inds)` — canonical relabeling of normal-mode indices; resonance
  conditions by `len(pf)`; properties by differentiation order then operator order.
- `VibPerturbedTerm.h(also_sort=True, nm_inds=...)`, consumed by `terms_simplify` to detect
  "algebraically identical up to coefficient".
- `PolProp.__eq__` **includes `self.inds`** — indices are part of identity.

**The evaluator's notion of identity**, same payload:

- `PropsCollection.sort()` — sorts by `ops[0].o`, the Cartesian axis, to push non-averaged factors last.
- `identify_avrg_motif()` — deep-copies and **sets `inds = None`**; indices deliberately erased.
- `PropsCollection.__hash__` = `(cart_axes, mode_indices)`.

Derive: *"the same term algebraically after canonical relabeling — indices matter."*
Evaluator: *"resolves to the same precomputed averaging tensor — indices are noise."*

One erases what the other keeps. **A class has one `__eq__`.** That is arithmetic, not preference.
`PropsCollection` therefore cannot be shared — not because derive wouldn't want a props container,
but because the container derive wants has a different sameness baked in.

Corroborating signal: `sort` denotes two unrelated operations at the two levels. When the same word
means different things in two layers, it is usually one fact surfacing as a naming collision.

---

## 4. "Wrapper" is not one abstraction — it is three

Calling them "wrappers for collections of symbolic abstractions" makes them sound unitary.
They are three families glued together, and the seam is the ownership line.

**(a) Predicates over the payload** — `get_cart_axes`, `get_mode_indices`, `get_total_difforder`,
and the `bool(p.ops)` test under the averaged/non-averaged split. Genuinely derive-level facts.
Most don't need the collection at all (`get_cart_axes` is a comprehension over `p.ops`).
If derive ever wants them: free functions over `Sequence[PolProp]`, or accessors on `PolProp`.

**(b) Identity and canonicalization** — `sort`, `__eq__`, `__hash__`, `_tuplify`,
`identify_avrg_motif`. Settled by §3. Evaluator's.

**(c) Lookup-key construction** — `get_num_indices_vibenedenom`, `get_mode_indices_group_template`,
`make_vibdiff_key`. These exist because of a chosen storage layout downstream. Evaluator's.

Precision on (a): *whether a factor carries Cartesian operators* is a physics predicate and is fine
in derive. What is not physics is **splitting the product into two disjoint sets returned as a pair** —
that shape exists to feed two lookup pathways (precomputed averaging tensor vs. direct array index).
The predicate goes down; the split stays up.

Separate the three and the ownership question dissolves: push (a) down as small pure accessors,
and what remains is a frozen key type with a canonical tuple, which nobody would argue about.

---

## 5. Mechanical rules

Not vibes — things that can be checked, ideally in CI.

### R1 — Derive types cross exactly one boundary

> One module in `wilson_intensities` imports `wilson_derive`. Below it, `VibPerturbedTerm`,
> `PolProp`, `VibDiffTerm`, `ResonanceCondition` do not appear.

Grep-able; enforceable with an import-linter contract or a three-line test.
Current count: **8 modules** in `amplitudes/` import `wilson_derive`. Target: 1.

### R2 — Identity ownership

Read `__eq__` / `__hash__` / the canonical form. If it is shaped by a downstream storage or dispatch
strategy, the type is downstream's. (§3.)

### R3 — Mutation direction

> A layer may not mutate objects owned by an upstream layer.

`_set_attr_for_all_props` assigns to `PolProp.inds` — a derive-owned object, possibly shared with
other holders of that term. `identify_avrg_motif` deep-copies first; `sort()` does not.
Grep-able (attribute assignment onto derive types inside intensities). This independently forces
translate-over-wrap: you cannot canonicalize by mutation if you do not own the object, and
canonicalizing by deep-copy is admitting you made a new type anyway.

### R4 — The strategy-change test

For each method ask: *if precomputed averaging tensors were replaced by on-the-fly averaging,
does this survive?* `identify_avrg_motif`, `sort`, `get_num_indices_vibenedenom` die.
`get_cart_axes` survives. A class whose rate of change tracks the evaluator's storage strategy
is the evaluator's class.

### R5 — The no-molecule test

The plan stage must run to completion with zero molecular data, and its output must be assertable.
If a stage cannot be tested without a molecule, it is not the plan stage.

### R6 — Asymmetric cost of speculative sharing

- Put them in derive, derive never uses them → derive carries API surface for a downstream consumer;
  the evaluator's canonicalization constrains derive's release surface; and when the evaluator needs
  a *different* equality, it is stuck (§3).
- Put them in intensities, derive later needs them → move the (a)-family accessors down. Small, pure,
  local, mechanical.

The second mistake is far cheaper. **Put shared-looking code at the consumer; promote on the second
real use.** Demoting is easy; retracting from a shared layer is not. (README rule 18.)

---

## 6. The missing stage: plan / compile

Contract:

- **in:** terms + axis choice. **No molecule.**
- **out:** `EvalPlan` — compiled terms plus a **work manifest**: which averaging motifs to build,
  which vibenedenom index sets to tabulate, which property trivial names are needed.
- **property:** pure function; cacheable, serializable, testable offline (R5).

Payoffs:

- `TermsInAxes.need_what()` (`rps_evaluation.py:65`, currently `NotImplementedError`) *is* the work
  manifest. It enables "verify all data present before computing anything", and batching of QC jobs
  from the manifest — which `MolPropsCollection.group_by_calc_setup` / `build_request_dict` are
  approaching from the wrong end.
- It removes existing duplication. `parse_vibpert_term` computes `avrg_expr`, `non_avrg_expr`,
  `freqterms`, `idx_summ/idx_nonsumm` — and `evaluate_term_coeffs` recomputes all four from the raw
  term, inside the loop, while `RspEvalTerm` is consumed by nothing. The parse exists but is unused.
  Make `RspEvalTerm` the only thing the evaluator sees and both problems close.

Vocabulary: three stages want three words. Derive's types are **symbolic**, plan's are **compiled**,
the numeric stage's are **evaluated**. Much of the confusion in the `# symbolic` comments comes from
having only two words for three things.

---

## 7. Wrap vs. translate

**Translate**, in the plan stage, into intensities-owned value types.

1. A canonical form is needed that derive does not provide — that is what `sort()` and index-stripping are.
2. `_set_attr_for_all_props` mutates derive-owned objects in place (R3). The deep-copy in
   `identify_avrg_motif` shows the hazard is already known; the cost of translation is already
   being paid, just defensively and inconsistently.
3. Derive can then evolve `PolProp` without silently breaking evaluation dict keys.
4. The derive boundary becomes literally one function, making R1 cheap to enforce.

Honest cost: one translation function to maintain, and a second vocabulary for the same physics.
If `PolProp` were stable *and* the evaluator needed no canonical form, wrapping would be less code.
Neither holds.

---

## 8. The coupling that matters more: `wilson_main`

`eval.py` lines ~355–650 contain verbatim copies of `DataOriginInfo`, `MolecularProperty`,
`MolPropsCollection`, `VibState` from `wilson_main/abstractions.py` — tab indentation included,
which is why the file is mixed tabs/spaces.

The copy is a symptom, not a fix. The situation:

- `wilson_main` is the orchestration layer and imports `amplitudes`, including at module level
  (`spectrum_abstractions.py:5`).
- `amplitudes` imports `wilson_main.abstractions` from 8 modules, some at module level
  (`evaluators.py:8`).

A package-level cycle that works only because the two modules that would close it happen not to meet.
Load-order fragile; the deferred function-body imports at `workflow_abstractions.py:397` are it
already biting.

Copying the classes into intensities makes intensities a *second home* for the molecular data model —
the opposite of the separation being sought. The fix is direction, not duplication: **the data model
is a leaf below both.** `MolecularProperty` / `MolPropsCollection` / `DataOriginInfo` / `VibState`
depend on nothing in the suite and should live in a package that both `wilson_main` and
`wilson_intensities` import and that imports neither. (`molprops.py` looks like an abandoned attempt
at this, except it holds *instances* and imports `wilson_main`, so it is on the wrong side.)

The one genuinely new class in that block: **`MolSystemData`** — "everything obtained externally,
no configuration." Best idea in the draft. Keep it, and make it the only door molecular data enters by.

Consequence: `EvaluationDataAndConfigs` mixes data (`props_data`, `vibstates_data`, `nc_sqrt_eigval`)
with configuration (`nm_inds_choices`, `pulse_polarization_vector`) — which is exactly why
`_prep_data` needs both `data` and `setup`. Split by who-reads-it (README rule 16): the plan stage
takes config only, the numeric stage takes data only. Both `EvaluationDataAndConfigs` and
`PrecalculatedData` are all-`None`-default grab-bags, which README rule 2 already forbids.

---

## 9. Value-type leaks in the draft

Same bug class the README's closing paragraph was written against.

**Asymmetric equality on dict keys.** `PropsCollection.__eq__` is
`all(p in other.props for p in self.props)` — ignores length and multiplicity, not symmetric.
`FreqTermsCollection.__eq__` has the same shape. These are keys in `avrg_expr_tensor_mapping`;
asymmetric `__eq__` on a hash key is a latent wrong-tensor lookup. `ResonanceMotif` gets it right
via `_tuplify()` — that is the pattern: one canonical tuple, both dunders derived from it.

**Mutation after key insertion.** `sort()` mutates `self.props`, returns `self`, and reassigns a
*list* into a field `__post_init__` had normalized to a tuple. A collection used as a key can be
mutated after insertion. Frozen, with `sorted()` returning a new instance.

**`ParameterSet` has a type that lies.** Declared `Mapping[str, int]`, but `__init__` injects
`params['zero'] = 'zero'` (a `str` value) and `__getitem__` remaps `''` → `'zero'`; `__lt__`
hardcodes the alphabet `('a'..'h')`. A generic index-assignment type that secretly knows
vibrational-state labelling conventions. Decide which it is: if generic, the zero sentinel and
ordering are policy living in a labelling module; if domain, name it (`IndexAssignment`) and make
the conventions explicit and tested. The ground state currently spelled three ways
(`''`, `'zero'`, `state_label == 'zero'`) is that ambiguity leaking.

**Caching threaded through domain types.** `VibDiff.cache_it(vibdiff_cache)` — a domain object
mutating a cache handed to it — plus `VibDiffCache` keyed on label strings with sign-flip-on-
normalized-miss. Violates README rules 5 and 7. The energy difference is a pure function of a
normalized label pair: make it one, memoize at function level, and `VibDiff` need not know a cache
exists. `make_vibdiff_key` then belongs to the plan stage, being a key function.

**Drafting slips worth noting** because they point at the unused-parse problem: `parse_vibpert_term`
never passes `term_id` (required field → `TypeError`), and `num_coeff` is typed
`FreqTermsCollection` but assigned a `float`.

---

## 10. Proposed shape

```
rsp_eval/
  plan.py        # ONLY importer of wilson_derive
                 #   in:  terms, axis choice        out: EvalPlan + WorkManifest
                 #   holds: PropsCollection, FreqTermsCollection, ResonanceMotif,
                 #          parse_vibpert_term, index bookkeeping, motif keys
  ingest.py      # MolSystemData, VibStatesData — the data door
  precompute.py  # manifest + data -> Tables (avrg tensors, vibenedenom, vibdiff energies)
  kernel.py      # arrays and floats only; the hierarchical sum; picklable
  features.py    # locations + coefficients -> SpectralFeature
  render.py      # grid
```

Arrows one-way, left to right. `plan.py` imports derive; nothing else does. `kernel.py` imports
nothing of ours above the array level (README rules 11 and 14) — which is what makes parallel
execution and testing-without-a-setup work.

Boundary types: `EvalPlan` (frozen, the compiled program), `WorkManifest` (what data/tensors are
needed; answers `need_what()`), `Tables` (numeric precomputation keyed by manifest entries),
`Coefficients` (`dict[IndexAssignment, complex]`), `SpectralFeature` (exists).

---

## 11. Open questions

1. **Is the plan a function of the axis choice, or only of the terms?** `TermsInAxes` bundles them,
   so changing axes invalidates the whole plan. If the property/index structure is axis-independent
   and only resonance location depends on axes, that is two sub-stages with different invalidation —
   which matters a lot when sweeping over axis choices.

2. **Where does motif handling split?** `ResonanceMotif` is a plan-level dedup key, but
   `process_resonance_motifs` needs `vibstates_data` — molecular. Motif *identification* is plan,
   motif *location* is numeric; one function currently does both. Splitting makes the
   "many terms, one motif" saving explicit rather than incidental.

3. **Does `wilson_main` keep owning the data model?** If yes, intensities depends on wilson_main and
   the cycle risk is accepted. If no, the model moves to a leaf and both depend on it. There is no
   third answer; the copy-paste in the draft is what trying to have one looks like.

---

## Through-line

Name the compile stage. Put every derive-facing type in it. Translate rather than wrap at that
boundary. Push the data model *down* rather than copying it *sideways*. Let identity — not
payload — decide ownership.
