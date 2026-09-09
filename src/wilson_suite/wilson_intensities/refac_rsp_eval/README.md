## Stages

1. Derivation — perturbation theory produces terms with free symbols. Depends on order and experiment type, not on any molecule. Produces magn_conds.
2. Axis translation — terms re-expressed in the chosen spectral axes. Still symbolic.
3. Feature computation — molecular data arrives; resonance positions and complex coefficients become numbers. Output is list[SpectralFeature].
4. Rendering — features become an array on a grid.

---

## Order of drafting

Write the public function signatures first — names, argument types, return types, nothing else. If you cannot write a signature without a placeholder type, that type is the actual open question and you should resolve it before writing bodies.

Write the boundary types second. These are the values passed between stages. They constrain everything downstream, so getting them wrong is expensive to undo.

Write the pipeline third as one flat script with no functions and no modules — every call in sequence, top to bottom, in a scratch file. Read it. If the linear version is hard to follow, splitting it into modules will hide the problem rather than fix it. Cut into functions only at points where the flat version has a natural break, and let the module boundaries follow the function boundaries.

## Value types

1. Every type has exactly one construction path. Conversions live in the producer. If a consumer needs to accept three input shapes, the producers are inconsistent — fix them. `_prep_terms` was this failure.

2. No optional field that the code dereferences unconditionally. `x: T | None = None` means the pipeline has a working path where `x` is absent. If it does not, the field is required.

3. Frozen everywhere, tuples not lists. This is not for safety in the abstract; it makes it impossible to attach a derived value to an object after construction, which is the mechanism of the bug you are designing against.

4. No derived value stored as a field. Two allowed forms: a `@property` when the computation is cheap and depends only on that object's own fields, or a separate return value passed forward when it is expensive. There is no third case.

5. No object holds a reference to a later-stage object. `EvaluatedResult` may hold the setup that produced it. The setup may not hold results.

6. Units in the field name, always: `gamma_cm1`, `gamma_au`, `energy_hartree`. The conversion is a property or a function, never an in-place mutation.

## Functions

7. Pure functions of their arguments. No reads of module-level state. Every value a function uses is a parameter or derived from one.

8. No helper that wraps a single call. If the wrapper adds no validation, no branching, and no transformation, inline it.

9. Type checks at the entry points only. Inside a stage, trust the types. `isinstance` in the interior means the boundary contract is not enforced.

10. Validate at construction. An invariant spanning two fields is checked in the `__post_init__` of the object holding both. If the fields live on different objects, they belong on the same one.

11. The numerical kernel takes arrays and floats. No config objects, no domain types. This is what makes it picklable for parallel execution and testable without building a setup.

## Modules

12. One reason to open a file. `grid.py` is opened when grid geometry is wrong; `features.py` when a coefficient is wrong. If a bug could plausibly send you to either, the split is in the wrong place.

13. Import direction is one-way and matches stage order. Draw it once and check it: `setup_types` ← `features`, `setup_types` ← `render`, `grid`/`kernel` ← `render`, and `features` never imports `render` or vice versa. Any cycle means two things you thought were separate are one thing.

14. `kernel.py` imports nothing from your own package above the array level.

## Parameters

15. A parameter belongs to the earliest stage that cannot produce its output without it. State the test as: if changing this value lets you reuse the previous stage's output unchanged, it belongs to the later stage.

16. Group parameters by the stage that reads them, not by topic. Four render parameters in one `RenderSetup` beats a general `Config` with fourteen fields where a reader cannot tell which four matter.

17. No parameter appears in two config objects. If it is needed in two places, it is passed forward, not duplicated.

## Minimality

18. Introduce an abstraction when there are two concrete instances of it, not before. One executor, one axis choice, one lineshape means write the concrete version.

19. No base classes, no registries, no plugin dispatch until you have a second implementation in hand.

20. Delete rather than comment out. Commented fields with live code referencing them was how `calc_config` broke.

## One thing to write down

State the invariant the design maintains, in prose, at the top of `setup_types.py`: no value is stored after being computed from another stored value; every derived quantity is either a property over immutable inputs or a return value. Any future change that fails that sentence reintroduces the stale-cache class of bug, and having it written means you can check a diff against it rather than remembering.