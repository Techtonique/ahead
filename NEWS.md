# version 0.6.0

- Include basic methods: mean forecast, median forecast, random walk forecast


# version 0.5.0

- add dropout regularization to `ridge2f`
- parallel execution for `type_pi == bootstrap` in `ridge2f` (done in R /!\, experimental)
- preallocate matrices for `type_forecast == recursive` in `ridge2f`


# version 0.4.2

- new attributes mean, lower bound, upper bound forecast as numpy arrays


# version 0.4.1

- use `get_frequency` to get series frequency as a number
- create a function `get_tscv_indices` for getting time series cross-validation indices
