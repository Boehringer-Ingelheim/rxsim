# Example 2: Two arm \| Interim \| Continuous \| t-test

``` r
library(rxsim)
```

This example extends [Example
1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md)
by introducing two named analyses — `interim` and `final` — that fire at
different enrollment milestones. Understanding how analysis names
propagate through
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
is the foundation for multi-analysis designs shown in later examples.

## Scenario

Capture scenario parameters. We will assume piece-wise linear
enrollment.

``` r
sample_size <- 100
arms        <- c("pbo", "trt")
allocation  <- c(1, 1)
delta       <- 0.2  # true treatment - placebo mean difference
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn    <- function(n) rexp(n, rate = 0.01)
scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation  = list(allocation),
  delta       = delta
)
```

`allocation = c(1, 1)` specifies balanced randomisation.
[`tidyr::expand_grid()`](https://tidyr.tidyverse.org/reference/expand_grid.html)
embeds design parameters directly into each result row for traceability
across parameter sweeps.

## Populations

Define population generators.

``` r
population_generators <- list(
  pbo = function(n) data.frame(
    id = 1:n,
    value = rnorm(n),
    readout_time = 1
  ),
  trt = function(n) data.frame(
    id = 1:n,
    value = rnorm(n, mean = delta),
    readout_time = 1
  )
)
```

Each generator is a function of `n` returning a `data.frame` with one
row per subject. `readout_time = 1` means the endpoint is observed 1
time unit after enrollment. The treatment arm has a mean shift of δ =
0.2 SD units over placebo.

## Triggers & Analysis

Final analysis at full enrollment.

``` r
analysis_generators <- list(
  interim = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= (!!sample_size) / 2
    ),
    analysis = function(df, timer){
      df_enrolled <- df |> subset(!is.na(enroll_time))
      tt <- t.test(value ~ arm, data = df_enrolled)
      data.frame(
        scenario,
        n_total = nrow(df_enrolled),
        mean_pbo = mean(df_enrolled$value[df_enrolled$arm == "pbo"]),
        mean_trt = mean(df_enrolled$value[df_enrolled$arm == "trt"]),
        p_value = unname(tt$p.value),
        stringsAsFactors = FALSE
      )
    }
  ),
  
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer){
      df_enrolled <- df |> subset(!is.na(enroll_time))
      tt <- t.test(value ~ arm, data = df_enrolled)
      data.frame(
        scenario,
        n_total = nrow(df_enrolled),
        mean_pbo = mean(df_enrolled$value[df_enrolled$arm == "pbo"]),
        mean_trt = mean(df_enrolled$value[df_enrolled$arm == "trt"]),
        p_value = unname(tt$p.value),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

Two analyses are defined. The `interim` analysis fires at half
enrollment and the `final` analysis fires when all subjects are
enrolled. These names appear in the `analysis` column of
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
output, making it straightforward to distinguish interim from final
results when stacking rows across timepoints.

## Trial

Make multiple trial replicates.

``` r
trials <- replicate_trial(
  trial_name = "test_trial",
  sample_size = sample_size,
  arms = arms,
  allocation = allocation,
  enrollment = enrollment_fn,
  dropout = dropout_fn,
  analysis_generators = analysis_generators,
  population_generators = population_generators,
  n = 3
)
```

## Simulate

To simulate all replicates:

``` r
run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: test_trial_1
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[2]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: test_trial_2
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[3]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: test_trial_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
```

Bind one row per replicate into one data frame.

``` r
replicate_results <- collect_results(trials)
replicate_results
#>     replicate  timepoint analysis sample_size allocation delta n_total
#> 1           1   48.42750  interim         100       1, 1   0.2      50
#> 2           1   49.85442  interim         100       1, 1   0.2      51
#> 3           1   49.95773  interim         100       1, 1   0.2      52
#> 4           1   51.72320  interim         100       1, 1   0.2      53
#> 5           1   55.51429  interim         100       1, 1   0.2      54
#> 6           1   56.14240  interim         100       1, 1   0.2      55
#> 7           1   57.45119  interim         100       1, 1   0.2      56
#> 8           1   59.12747  interim         100       1, 1   0.2      57
#> 9           1   59.14277  interim         100       1, 1   0.2      58
#> 10          1   59.46349  interim         100       1, 1   0.2      59
#> 11          1   59.48707  interim         100       1, 1   0.2      60
#> 12          1   60.15817  interim         100       1, 1   0.2      61
#> 13          1   60.57574  interim         100       1, 1   0.2      62
#> 14          1   60.58369  interim         100       1, 1   0.2      63
#> 15          1   60.83251  interim         100       1, 1   0.2      64
#> 16          1   61.10046  interim         100       1, 1   0.2      65
#> 17          1   61.27583  interim         100       1, 1   0.2      66
#> 18          1   61.44810  interim         100       1, 1   0.2      67
#> 19          1   62.75277  interim         100       1, 1   0.2      68
#> 20          1   63.29203  interim         100       1, 1   0.2      69
#> 21          1   65.48555  interim         100       1, 1   0.2      70
#> 22          1   65.81196  interim         100       1, 1   0.2      71
#> 23          1   66.19484  interim         100       1, 1   0.2      72
#> 24          1   66.34897  interim         100       1, 1   0.2      73
#> 25          1   66.72386  interim         100       1, 1   0.2      74
#> 26          1   68.14108  interim         100       1, 1   0.2      75
#> 27          1   68.16704  interim         100       1, 1   0.2      76
#> 28          1   68.35415  interim         100       1, 1   0.2      77
#> 29          1   70.30805  interim         100       1, 1   0.2      78
#> 30          1   70.70770  interim         100       1, 1   0.2      79
#> 31          1   72.24447  interim         100       1, 1   0.2      80
#> 32          1   73.16830  interim         100       1, 1   0.2      81
#> 33          1   73.87027  interim         100       1, 1   0.2      82
#> 34          1   75.19949  interim         100       1, 1   0.2      83
#> 35          1   75.83765  interim         100       1, 1   0.2      84
#> 36          1   76.24062  interim         100       1, 1   0.2      85
#> 37          1   76.38961  interim         100       1, 1   0.2      86
#> 38          1   76.39273  interim         100       1, 1   0.2      87
#> 39          1   76.69642  interim         100       1, 1   0.2      88
#> 40          1   77.03205  interim         100       1, 1   0.2      89
#> 41          1   78.35709  interim         100       1, 1   0.2      90
#> 42          1   79.45524  interim         100       1, 1   0.2      91
#> 43          1   79.93142  interim         100       1, 1   0.2      92
#> 44          1   80.25999  interim         100       1, 1   0.2      93
#> 45          1   82.70303  interim         100       1, 1   0.2      94
#> 46          1   82.75639  interim         100       1, 1   0.2      95
#> 47          1   85.00719  interim         100       1, 1   0.2      96
#> 48          1   85.38033  interim         100       1, 1   0.2      97
#> 49          1   87.07561  interim         100       1, 1   0.2      98
#> 50          1   87.29870  interim         100       1, 1   0.2      99
#> 51          1   87.36847  interim         100       1, 1   0.2     100
#> 52          1   87.36847    final         100       1, 1   0.2     100
#> 53          1  347.10044  interim         100       1, 1   0.2     100
#> 54          1  347.10044    final         100       1, 1   0.2     100
#> 55          1  364.07755  interim         100       1, 1   0.2     100
#> 56          1  364.07755    final         100       1, 1   0.2     100
#> 57          1  448.68714  interim         100       1, 1   0.2     100
#> 58          1  448.68714    final         100       1, 1   0.2     100
#> 59          1  469.18697  interim         100       1, 1   0.2     100
#> 60          1  469.18697    final         100       1, 1   0.2     100
#> 61          1  487.15448  interim         100       1, 1   0.2     100
#> 62          1  487.15448    final         100       1, 1   0.2     100
#> 63          1  538.98616  interim         100       1, 1   0.2     100
#> 64          1  538.98616    final         100       1, 1   0.2     100
#> 65          1  606.20122  interim         100       1, 1   0.2     100
#> 66          1  606.20122    final         100       1, 1   0.2     100
#> 67          1  658.76511  interim         100       1, 1   0.2     100
#> 68          1  658.76511    final         100       1, 1   0.2     100
#> 69          1  794.98780  interim         100       1, 1   0.2     100
#> 70          1  794.98780    final         100       1, 1   0.2     100
#> 71          1  944.07711  interim         100       1, 1   0.2     100
#> 72          1  944.07711    final         100       1, 1   0.2     100
#> 73          1 1181.37320  interim         100       1, 1   0.2     100
#> 74          1 1181.37320    final         100       1, 1   0.2     100
#> 75          1 1212.56973  interim         100       1, 1   0.2     100
#> 76          1 1212.56973    final         100       1, 1   0.2     100
#> 77          1 1232.97050  interim         100       1, 1   0.2     100
#> 78          1 1232.97050    final         100       1, 1   0.2     100
#> 79          1 1264.36967  interim         100       1, 1   0.2     100
#> 80          1 1264.36967    final         100       1, 1   0.2     100
#> 81          1 1365.41125  interim         100       1, 1   0.2     100
#> 82          1 1365.41125    final         100       1, 1   0.2     100
#> 83          1 1414.98418  interim         100       1, 1   0.2     100
#> 84          1 1414.98418    final         100       1, 1   0.2     100
#> 85          1 1580.61613  interim         100       1, 1   0.2     100
#> 86          1 1580.61613    final         100       1, 1   0.2     100
#> 87          1 1905.27865  interim         100       1, 1   0.2     100
#> 88          1 1905.27865    final         100       1, 1   0.2     100
#> 89          1 1928.66156  interim         100       1, 1   0.2     100
#> 90          1 1928.66156    final         100       1, 1   0.2     100
#> 91          1 1948.35656  interim         100       1, 1   0.2     100
#> 92          1 1948.35656    final         100       1, 1   0.2     100
#> 93          1 2080.41273  interim         100       1, 1   0.2     100
#> 94          1 2080.41273    final         100       1, 1   0.2     100
#> 95          1 2116.29566  interim         100       1, 1   0.2     100
#> 96          1 2116.29566    final         100       1, 1   0.2     100
#> 97          1 2471.76965  interim         100       1, 1   0.2     100
#> 98          1 2471.76965    final         100       1, 1   0.2     100
#> 99          1 2506.45921  interim         100       1, 1   0.2     100
#> 100         1 2506.45921    final         100       1, 1   0.2     100
#> 101         1 2624.28188  interim         100       1, 1   0.2     100
#> 102         1 2624.28188    final         100       1, 1   0.2     100
#> 103         1 2670.38674  interim         100       1, 1   0.2     100
#> 104         1 2670.38674    final         100       1, 1   0.2     100
#> 105         1 2721.59496  interim         100       1, 1   0.2     100
#> 106         1 2721.59496    final         100       1, 1   0.2     100
#> 107         1 2730.16239  interim         100       1, 1   0.2     100
#> 108         1 2730.16239    final         100       1, 1   0.2     100
#> 109         1 2897.71439  interim         100       1, 1   0.2     100
#> 110         1 2897.71439    final         100       1, 1   0.2     100
#> 111         1 2943.10416  interim         100       1, 1   0.2     100
#> 112         1 2943.10416    final         100       1, 1   0.2     100
#> 113         1 2970.64127  interim         100       1, 1   0.2     100
#> 114         1 2970.64127    final         100       1, 1   0.2     100
#> 115         1 3098.51997  interim         100       1, 1   0.2     100
#> 116         1 3098.51997    final         100       1, 1   0.2     100
#> 117         1 3101.99965  interim         100       1, 1   0.2     100
#> 118         1 3101.99965    final         100       1, 1   0.2     100
#> 119         1 3117.23713  interim         100       1, 1   0.2     100
#> 120         1 3117.23713    final         100       1, 1   0.2     100
#> 121         1 3264.59026  interim         100       1, 1   0.2     100
#> 122         1 3264.59026    final         100       1, 1   0.2     100
#> 123         1 3324.84528  interim         100       1, 1   0.2     100
#> 124         1 3324.84528    final         100       1, 1   0.2     100
#> 125         1 3376.45601  interim         100       1, 1   0.2     100
#> 126         1 3376.45601    final         100       1, 1   0.2     100
#> 127         1 3382.96904  interim         100       1, 1   0.2     100
#> 128         1 3382.96904    final         100       1, 1   0.2     100
#> 129         1 3392.33000  interim         100       1, 1   0.2     100
#> 130         1 3392.33000    final         100       1, 1   0.2     100
#> 131         1 3653.75656  interim         100       1, 1   0.2     100
#> 132         1 3653.75656    final         100       1, 1   0.2     100
#> 133         1 3778.41117  interim         100       1, 1   0.2     100
#> 134         1 3778.41117    final         100       1, 1   0.2     100
#> 135         1 3954.92212  interim         100       1, 1   0.2     100
#> 136         1 3954.92212    final         100       1, 1   0.2     100
#> 137         1 3993.06729  interim         100       1, 1   0.2     100
#> 138         1 3993.06729    final         100       1, 1   0.2     100
#> 139         1 4028.10899  interim         100       1, 1   0.2     100
#> 140         1 4028.10899    final         100       1, 1   0.2     100
#> 141         1 4041.70989  interim         100       1, 1   0.2     100
#> 142         1 4041.70989    final         100       1, 1   0.2     100
#> 143         1 4165.86312  interim         100       1, 1   0.2     100
#> 144         1 4165.86312    final         100       1, 1   0.2     100
#> 145         1 4195.93190  interim         100       1, 1   0.2     100
#> 146         1 4195.93190    final         100       1, 1   0.2     100
#> 147         1 4258.85114  interim         100       1, 1   0.2     100
#> 148         1 4258.85114    final         100       1, 1   0.2     100
#> 149         1 4480.33893  interim         100       1, 1   0.2     100
#> 150         1 4480.33893    final         100       1, 1   0.2     100
#> 151         1 4485.70500  interim         100       1, 1   0.2     100
#> 152         1 4485.70500    final         100       1, 1   0.2     100
#> 153         1 4538.39996  interim         100       1, 1   0.2     100
#> 154         1 4538.39996    final         100       1, 1   0.2     100
#> 155         1 4621.85701  interim         100       1, 1   0.2     100
#> 156         1 4621.85701    final         100       1, 1   0.2     100
#> 157         1 4980.81111  interim         100       1, 1   0.2     100
#> 158         1 4980.81111    final         100       1, 1   0.2     100
#> 159         1 5022.77206  interim         100       1, 1   0.2     100
#> 160         1 5022.77206    final         100       1, 1   0.2     100
#> 161         1 5025.26114  interim         100       1, 1   0.2     100
#> 162         1 5025.26114    final         100       1, 1   0.2     100
#> 163         1 5283.65388  interim         100       1, 1   0.2     100
#> 164         1 5283.65388    final         100       1, 1   0.2     100
#> 165         1 5322.97061  interim         100       1, 1   0.2     100
#> 166         1 5322.97061    final         100       1, 1   0.2     100
#> 167         1 5454.57682  interim         100       1, 1   0.2     100
#> 168         1 5454.57682    final         100       1, 1   0.2     100
#> 169         1 5667.53959  interim         100       1, 1   0.2     100
#> 170         1 5667.53959    final         100       1, 1   0.2     100
#> 171         1 5807.36313  interim         100       1, 1   0.2     100
#> 172         1 5807.36313    final         100       1, 1   0.2     100
#> 173         1 5822.36146  interim         100       1, 1   0.2     100
#> 174         1 5822.36146    final         100       1, 1   0.2     100
#> 175         1 6024.11144  interim         100       1, 1   0.2     100
#> 176         1 6024.11144    final         100       1, 1   0.2     100
#> 177         1 6033.82973  interim         100       1, 1   0.2     100
#> 178         1 6033.82973    final         100       1, 1   0.2     100
#> 179         1 6099.33001  interim         100       1, 1   0.2     100
#> 180         1 6099.33001    final         100       1, 1   0.2     100
#> 181         1 6240.00557  interim         100       1, 1   0.2     100
#> 182         1 6240.00557    final         100       1, 1   0.2     100
#> 183         1 6313.79473  interim         100       1, 1   0.2     100
#> 184         1 6313.79473    final         100       1, 1   0.2     100
#> 185         1 6398.47626  interim         100       1, 1   0.2     100
#> 186         1 6398.47626    final         100       1, 1   0.2     100
#> 187         1 6843.06418  interim         100       1, 1   0.2     100
#> 188         1 6843.06418    final         100       1, 1   0.2     100
#> 189         1 6893.26695  interim         100       1, 1   0.2     100
#> 190         1 6893.26695    final         100       1, 1   0.2     100
#> 191         1 6937.06974  interim         100       1, 1   0.2     100
#> 192         1 6937.06974    final         100       1, 1   0.2     100
#> 193         1 6940.14080  interim         100       1, 1   0.2     100
#> 194         1 6940.14080    final         100       1, 1   0.2     100
#> 195         1 6951.95172  interim         100       1, 1   0.2     100
#> 196         1 6951.95172    final         100       1, 1   0.2     100
#> 197         1 7300.51551  interim         100       1, 1   0.2     100
#> 198         1 7300.51551    final         100       1, 1   0.2     100
#> 199         1 7357.70574  interim         100       1, 1   0.2     100
#> 200         1 7357.70574    final         100       1, 1   0.2     100
#> 201         1 7422.76732  interim         100       1, 1   0.2     100
#> 202         1 7422.76732    final         100       1, 1   0.2     100
#> 203         1 7434.85486  interim         100       1, 1   0.2     100
#> 204         1 7434.85486    final         100       1, 1   0.2     100
#> 205         1 7455.60494  interim         100       1, 1   0.2     100
#> 206         1 7455.60494    final         100       1, 1   0.2     100
#> 207         1 7516.69768  interim         100       1, 1   0.2     100
#> 208         1 7516.69768    final         100       1, 1   0.2     100
#> 209         1 7811.33067  interim         100       1, 1   0.2     100
#> 210         1 7811.33067    final         100       1, 1   0.2     100
#> 211         1 7857.91280  interim         100       1, 1   0.2     100
#> 212         1 7857.91280    final         100       1, 1   0.2     100
#> 213         1 7997.65149  interim         100       1, 1   0.2     100
#> 214         1 7997.65149    final         100       1, 1   0.2     100
#> 215         1 8034.86411  interim         100       1, 1   0.2     100
#> 216         1 8034.86411    final         100       1, 1   0.2     100
#> 217         1 8063.25236  interim         100       1, 1   0.2     100
#> 218         1 8063.25236    final         100       1, 1   0.2     100
#> 219         1 8163.50249  interim         100       1, 1   0.2     100
#> 220         1 8163.50249    final         100       1, 1   0.2     100
#> 221         1 8287.81000  interim         100       1, 1   0.2     100
#> 222         1 8287.81000    final         100       1, 1   0.2     100
#> 223         1 8329.31008  interim         100       1, 1   0.2     100
#> 224         1 8329.31008    final         100       1, 1   0.2     100
#> 225         1 8368.84116  interim         100       1, 1   0.2     100
#> 226         1 8368.84116    final         100       1, 1   0.2     100
#> 227         1 8400.65785  interim         100       1, 1   0.2     100
#> 228         1 8400.65785    final         100       1, 1   0.2     100
#> 229         1 8571.29608  interim         100       1, 1   0.2     100
#> 230         1 8571.29608    final         100       1, 1   0.2     100
#> 231         1 8584.67803  interim         100       1, 1   0.2     100
#> 232         1 8584.67803    final         100       1, 1   0.2     100
#> 233         1 8625.88371  interim         100       1, 1   0.2     100
#> 234         1 8625.88371    final         100       1, 1   0.2     100
#> 235         1 8644.56538  interim         100       1, 1   0.2     100
#> 236         1 8644.56538    final         100       1, 1   0.2     100
#> 237         1 8700.35951  interim         100       1, 1   0.2     100
#> 238         1 8700.35951    final         100       1, 1   0.2     100
#> 239         1 8828.78289  interim         100       1, 1   0.2     100
#> 240         1 8828.78289    final         100       1, 1   0.2     100
#> 241         1 8850.81043  interim         100       1, 1   0.2     100
#> 242         1 8850.81043    final         100       1, 1   0.2     100
#> 243         1 9065.76198  interim         100       1, 1   0.2     100
#> 244         1 9065.76198    final         100       1, 1   0.2     100
#> 245         1 9187.05539  interim         100       1, 1   0.2     100
#> 246         1 9187.05539    final         100       1, 1   0.2     100
#> 247         2   40.26049  interim         100       1, 1   0.2      50
#> 248         2   40.41300  interim         100       1, 1   0.2      51
#> 249         2   42.80177  interim         100       1, 1   0.2      52
#> 250         2   43.97187  interim         100       1, 1   0.2      53
#> 251         2   45.37149  interim         100       1, 1   0.2      54
#> 252         2   45.85124  interim         100       1, 1   0.2      55
#> 253         2   46.17569  interim         100       1, 1   0.2      56
#> 254         2   46.25208  interim         100       1, 1   0.2      57
#> 255         2   47.26904  interim         100       1, 1   0.2      58
#> 256         2   47.52640  interim         100       1, 1   0.2      59
#> 257         2   48.76312  interim         100       1, 1   0.2      60
#> 258         2   50.25389  interim         100       1, 1   0.2      61
#> 259         2   51.07229  interim         100       1, 1   0.2      62
#> 260         2   51.48618  interim         100       1, 1   0.2      63
#> 261         2   52.64907  interim         100       1, 1   0.2      64
#> 262         2   54.82302  interim         100       1, 1   0.2      65
#> 263         2   54.99248  interim         100       1, 1   0.2      66
#> 264         2   54.99797  interim         100       1, 1   0.2      67
#> 265         2   55.24091  interim         100       1, 1   0.2      68
#> 266         2   56.14004  interim         100       1, 1   0.2      69
#> 267         2   56.99578  interim         100       1, 1   0.2      70
#> 268         2   57.12977  interim         100       1, 1   0.2      71
#> 269         2   58.28487  interim         100       1, 1   0.2      72
#> 270         2   58.98881  interim         100       1, 1   0.2      73
#> 271         2   59.07894  interim         100       1, 1   0.2      74
#> 272         2   60.19569  interim         100       1, 1   0.2      75
#> 273         2   60.27796  interim         100       1, 1   0.2      76
#> 274         2   60.32526  interim         100       1, 1   0.2      77
#> 275         2   60.99938  interim         100       1, 1   0.2      78
#> 276         2   63.84410  interim         100       1, 1   0.2      79
#> 277         2   67.00320  interim         100       1, 1   0.2      80
#> 278         2   67.58740  interim         100       1, 1   0.2      81
#> 279         2   69.59529  interim         100       1, 1   0.2      82
#> 280         2   71.77546  interim         100       1, 1   0.2      83
#> 281         2   71.82889  interim         100       1, 1   0.2      84
#> 282         2   71.99695  interim         100       1, 1   0.2      85
#> 283         2   72.25126  interim         100       1, 1   0.2      86
#> 284         2   73.98488  interim         100       1, 1   0.2      87
#> 285         2   74.76001  interim         100       1, 1   0.2      88
#> 286         2   75.92301  interim         100       1, 1   0.2      89
#> 287         2   76.10327  interim         100       1, 1   0.2      90
#> 288         2   76.41787  interim         100       1, 1   0.2      91
#> 289         2   78.71832  interim         100       1, 1   0.2      92
#> 290         2   79.13490  interim         100       1, 1   0.2      93
#> 291         2   79.20720  interim         100       1, 1   0.2      94
#> 292         2   79.89032  interim         100       1, 1   0.2      95
#> 293         2   81.04052  interim         100       1, 1   0.2      96
#> 294         2   81.28636  interim         100       1, 1   0.2      97
#> 295         2   82.09278  interim         100       1, 1   0.2      98
#> 296         2   83.62521  interim         100       1, 1   0.2      99
#> 297         2   83.92580  interim         100       1, 1   0.2     100
#> 298         2   83.92580    final         100       1, 1   0.2     100
#> 299         2  162.86506  interim         100       1, 1   0.2     100
#> 300         2  162.86506    final         100       1, 1   0.2     100
#> 301         2  164.52703  interim         100       1, 1   0.2     100
#> 302         2  164.52703    final         100       1, 1   0.2     100
#> 303         2  409.86390  interim         100       1, 1   0.2     100
#> 304         2  409.86390    final         100       1, 1   0.2     100
#> 305         2  410.46923  interim         100       1, 1   0.2     100
#> 306         2  410.46923    final         100       1, 1   0.2     100
#> 307         2  711.15813  interim         100       1, 1   0.2     100
#> 308         2  711.15813    final         100       1, 1   0.2     100
#> 309         2  775.96471  interim         100       1, 1   0.2     100
#> 310         2  775.96471    final         100       1, 1   0.2     100
#> 311         2  859.67365  interim         100       1, 1   0.2     100
#> 312         2  859.67365    final         100       1, 1   0.2     100
#> 313         2 1025.88990  interim         100       1, 1   0.2     100
#> 314         2 1025.88990    final         100       1, 1   0.2     100
#> 315         2 1131.25232  interim         100       1, 1   0.2     100
#> 316         2 1131.25232    final         100       1, 1   0.2     100
#> 317         2 1156.72552  interim         100       1, 1   0.2     100
#> 318         2 1156.72552    final         100       1, 1   0.2     100
#> 319         2 1221.70859  interim         100       1, 1   0.2     100
#> 320         2 1221.70859    final         100       1, 1   0.2     100
#> 321         2 1376.43033  interim         100       1, 1   0.2     100
#> 322         2 1376.43033    final         100       1, 1   0.2     100
#> 323         2 1382.39410  interim         100       1, 1   0.2     100
#> 324         2 1382.39410    final         100       1, 1   0.2     100
#> 325         2 1383.07312  interim         100       1, 1   0.2     100
#> 326         2 1383.07312    final         100       1, 1   0.2     100
#> 327         2 1452.59714  interim         100       1, 1   0.2     100
#> 328         2 1452.59714    final         100       1, 1   0.2     100
#> 329         2 1458.82499  interim         100       1, 1   0.2     100
#> 330         2 1458.82499    final         100       1, 1   0.2     100
#> 331         2 1487.90153  interim         100       1, 1   0.2     100
#> 332         2 1487.90153    final         100       1, 1   0.2     100
#> 333         2 1574.45724  interim         100       1, 1   0.2     100
#> 334         2 1574.45724    final         100       1, 1   0.2     100
#> 335         2 1687.75601  interim         100       1, 1   0.2     100
#> 336         2 1687.75601    final         100       1, 1   0.2     100
#> 337         2 1885.57444  interim         100       1, 1   0.2     100
#> 338         2 1885.57444    final         100       1, 1   0.2     100
#> 339         2 1930.83150  interim         100       1, 1   0.2     100
#> 340         2 1930.83150    final         100       1, 1   0.2     100
#> 341         2 2103.96690  interim         100       1, 1   0.2     100
#> 342         2 2103.96690    final         100       1, 1   0.2     100
#> 343         2 2266.56908  interim         100       1, 1   0.2     100
#> 344         2 2266.56908    final         100       1, 1   0.2     100
#> 345         2 2482.22496  interim         100       1, 1   0.2     100
#> 346         2 2482.22496    final         100       1, 1   0.2     100
#> 347         2 2829.74784  interim         100       1, 1   0.2     100
#> 348         2 2829.74784    final         100       1, 1   0.2     100
#> 349         2 2876.89371  interim         100       1, 1   0.2     100
#> 350         2 2876.89371    final         100       1, 1   0.2     100
#> 351         2 2904.98490  interim         100       1, 1   0.2     100
#> 352         2 2904.98490    final         100       1, 1   0.2     100
#> 353         2 2982.17747  interim         100       1, 1   0.2     100
#> 354         2 2982.17747    final         100       1, 1   0.2     100
#> 355         2 3051.66387  interim         100       1, 1   0.2     100
#> 356         2 3051.66387    final         100       1, 1   0.2     100
#> 357         2 3207.73134  interim         100       1, 1   0.2     100
#> 358         2 3207.73134    final         100       1, 1   0.2     100
#> 359         2 3233.42137  interim         100       1, 1   0.2     100
#> 360         2 3233.42137    final         100       1, 1   0.2     100
#> 361         2 3329.18229  interim         100       1, 1   0.2     100
#> 362         2 3329.18229    final         100       1, 1   0.2     100
#> 363         2 3335.55715  interim         100       1, 1   0.2     100
#> 364         2 3335.55715    final         100       1, 1   0.2     100
#> 365         2 3383.11618  interim         100       1, 1   0.2     100
#> 366         2 3383.11618    final         100       1, 1   0.2     100
#> 367         2 3397.06429  interim         100       1, 1   0.2     100
#> 368         2 3397.06429    final         100       1, 1   0.2     100
#> 369         2 3436.36201  interim         100       1, 1   0.2     100
#> 370         2 3436.36201    final         100       1, 1   0.2     100
#> 371         2 3476.51551  interim         100       1, 1   0.2     100
#> 372         2 3476.51551    final         100       1, 1   0.2     100
#> 373         2 3781.58448  interim         100       1, 1   0.2     100
#> 374         2 3781.58448    final         100       1, 1   0.2     100
#> 375         2 3890.17196  interim         100       1, 1   0.2     100
#> 376         2 3890.17196    final         100       1, 1   0.2     100
#> 377         2 3907.82140  interim         100       1, 1   0.2     100
#> 378         2 3907.82140    final         100       1, 1   0.2     100
#> 379         2 3957.68072  interim         100       1, 1   0.2     100
#> 380         2 3957.68072    final         100       1, 1   0.2     100
#> 381         2 4046.83243  interim         100       1, 1   0.2     100
#> 382         2 4046.83243    final         100       1, 1   0.2     100
#> 383         2 4070.13216  interim         100       1, 1   0.2     100
#> 384         2 4070.13216    final         100       1, 1   0.2     100
#> 385         2 4236.18139  interim         100       1, 1   0.2     100
#> 386         2 4236.18139    final         100       1, 1   0.2     100
#> 387         2 4366.13531  interim         100       1, 1   0.2     100
#> 388         2 4366.13531    final         100       1, 1   0.2     100
#> 389         2 5058.56822  interim         100       1, 1   0.2     100
#> 390         2 5058.56822    final         100       1, 1   0.2     100
#> 391         2 5195.06935  interim         100       1, 1   0.2     100
#> 392         2 5195.06935    final         100       1, 1   0.2     100
#> 393         2 5247.85622  interim         100       1, 1   0.2     100
#> 394         2 5247.85622    final         100       1, 1   0.2     100
#> 395         2 5306.75175  interim         100       1, 1   0.2     100
#> 396         2 5306.75175    final         100       1, 1   0.2     100
#> 397         2 5447.39854  interim         100       1, 1   0.2     100
#> 398         2 5447.39854    final         100       1, 1   0.2     100
#> 399         2 5453.48043  interim         100       1, 1   0.2     100
#> 400         2 5453.48043    final         100       1, 1   0.2     100
#> 401         2 5578.89134  interim         100       1, 1   0.2     100
#> 402         2 5578.89134    final         100       1, 1   0.2     100
#> 403         2 5613.89604  interim         100       1, 1   0.2     100
#> 404         2 5613.89604    final         100       1, 1   0.2     100
#> 405         2 5872.18930  interim         100       1, 1   0.2     100
#> 406         2 5872.18930    final         100       1, 1   0.2     100
#> 407         2 6068.78939  interim         100       1, 1   0.2     100
#> 408         2 6068.78939    final         100       1, 1   0.2     100
#> 409         2 6175.40092  interim         100       1, 1   0.2     100
#> 410         2 6175.40092    final         100       1, 1   0.2     100
#> 411         2 6190.98047  interim         100       1, 1   0.2     100
#> 412         2 6190.98047    final         100       1, 1   0.2     100
#> 413         2 6229.16404  interim         100       1, 1   0.2     100
#> 414         2 6229.16404    final         100       1, 1   0.2     100
#> 415         2 6323.41052  interim         100       1, 1   0.2     100
#> 416         2 6323.41052    final         100       1, 1   0.2     100
#> 417         2 6401.16754  interim         100       1, 1   0.2     100
#> 418         2 6401.16754    final         100       1, 1   0.2     100
#> 419         2 6496.77482  interim         100       1, 1   0.2     100
#> 420         2 6496.77482    final         100       1, 1   0.2     100
#> 421         2 6560.01235  interim         100       1, 1   0.2     100
#> 422         2 6560.01235    final         100       1, 1   0.2     100
#> 423         2 6569.19544  interim         100       1, 1   0.2     100
#> 424         2 6569.19544    final         100       1, 1   0.2     100
#> 425         2 6571.86589  interim         100       1, 1   0.2     100
#> 426         2 6571.86589    final         100       1, 1   0.2     100
#> 427         2 6807.20046  interim         100       1, 1   0.2     100
#> 428         2 6807.20046    final         100       1, 1   0.2     100
#> 429         2 6896.87699  interim         100       1, 1   0.2     100
#> 430         2 6896.87699    final         100       1, 1   0.2     100
#> 431         2 7109.88073  interim         100       1, 1   0.2     100
#> 432         2 7109.88073    final         100       1, 1   0.2     100
#> 433         2 7124.19290  interim         100       1, 1   0.2     100
#> 434         2 7124.19290    final         100       1, 1   0.2     100
#> 435         2 7159.79051  interim         100       1, 1   0.2     100
#> 436         2 7159.79051    final         100       1, 1   0.2     100
#> 437         2 7197.93433  interim         100       1, 1   0.2     100
#> 438         2 7197.93433    final         100       1, 1   0.2     100
#> 439         2 7316.61338  interim         100       1, 1   0.2     100
#> 440         2 7316.61338    final         100       1, 1   0.2     100
#> 441         2 7581.32242  interim         100       1, 1   0.2     100
#> 442         2 7581.32242    final         100       1, 1   0.2     100
#> 443         2 7847.68211  interim         100       1, 1   0.2     100
#> 444         2 7847.68211    final         100       1, 1   0.2     100
#> 445         2 7902.64242  interim         100       1, 1   0.2     100
#> 446         2 7902.64242    final         100       1, 1   0.2     100
#> 447         2 8062.46831  interim         100       1, 1   0.2     100
#> 448         2 8062.46831    final         100       1, 1   0.2     100
#> 449         2 8189.62566  interim         100       1, 1   0.2     100
#> 450         2 8189.62566    final         100       1, 1   0.2     100
#> 451         2 8233.68681  interim         100       1, 1   0.2     100
#> 452         2 8233.68681    final         100       1, 1   0.2     100
#> 453         2 8338.28947  interim         100       1, 1   0.2     100
#> 454         2 8338.28947    final         100       1, 1   0.2     100
#> 455         2 8346.03816  interim         100       1, 1   0.2     100
#> 456         2 8346.03816    final         100       1, 1   0.2     100
#> 457         2 8349.95850  interim         100       1, 1   0.2     100
#> 458         2 8349.95850    final         100       1, 1   0.2     100
#> 459         2 8359.23202  interim         100       1, 1   0.2     100
#> 460         2 8359.23202    final         100       1, 1   0.2     100
#> 461         2 8429.83124  interim         100       1, 1   0.2     100
#> 462         2 8429.83124    final         100       1, 1   0.2     100
#> 463         2 8441.34870  interim         100       1, 1   0.2     100
#> 464         2 8441.34870    final         100       1, 1   0.2     100
#> 465         2 8555.06586  interim         100       1, 1   0.2     100
#> 466         2 8555.06586    final         100       1, 1   0.2     100
#> 467         2 8585.46590  interim         100       1, 1   0.2     100
#> 468         2 8585.46590    final         100       1, 1   0.2     100
#> 469         2 8614.05634  interim         100       1, 1   0.2     100
#> 470         2 8614.05634    final         100       1, 1   0.2     100
#> 471         2 8780.06920  interim         100       1, 1   0.2     100
#> 472         2 8780.06920    final         100       1, 1   0.2     100
#> 473         2 8849.20991  interim         100       1, 1   0.2     100
#> 474         2 8849.20991    final         100       1, 1   0.2     100
#> 475         2 8851.85363  interim         100       1, 1   0.2     100
#> 476         2 8851.85363    final         100       1, 1   0.2     100
#> 477         2 8867.06403  interim         100       1, 1   0.2     100
#> 478         2 8867.06403    final         100       1, 1   0.2     100
#> 479         2 8968.75154  interim         100       1, 1   0.2     100
#> 480         2 8968.75154    final         100       1, 1   0.2     100
#> 481         2 9149.58844  interim         100       1, 1   0.2     100
#> 482         2 9149.58844    final         100       1, 1   0.2     100
#> 483         2 9263.90968  interim         100       1, 1   0.2     100
#> 484         2 9263.90968    final         100       1, 1   0.2     100
#> 485         2 9389.80570  interim         100       1, 1   0.2     100
#> 486         2 9389.80570    final         100       1, 1   0.2     100
#> 487         2 9399.46535  interim         100       1, 1   0.2     100
#> 488         2 9399.46535    final         100       1, 1   0.2     100
#> 489         2 9476.42755  interim         100       1, 1   0.2     100
#> 490         2 9476.42755    final         100       1, 1   0.2     100
#> 491         2 9483.19609  interim         100       1, 1   0.2     100
#> 492         2 9483.19609    final         100       1, 1   0.2     100
#> 493         2 9531.04625  interim         100       1, 1   0.2     100
#> 494         2 9531.04625    final         100       1, 1   0.2     100
#> 495         2 9574.95081  interim         100       1, 1   0.2     100
#> 496         2 9574.95081    final         100       1, 1   0.2     100
#> 497         2 9650.05829  interim         100       1, 1   0.2     100
#> 498         2 9650.05829    final         100       1, 1   0.2     100
#> 499         3   43.51109  interim         100       1, 1   0.2      50
#> 500         3   45.27401  interim         100       1, 1   0.2      51
#> 501         3   45.49279  interim         100       1, 1   0.2      52
#> 502         3   46.11225  interim         100       1, 1   0.2      53
#> 503         3   46.19910  interim         100       1, 1   0.2      54
#> 504         3   47.60653  interim         100       1, 1   0.2      55
#> 505         3   47.86172  interim         100       1, 1   0.2      56
#> 506         3   50.14635  interim         100       1, 1   0.2      57
#> 507         3   50.98130  interim         100       1, 1   0.2      58
#> 508         3   51.91686  interim         100       1, 1   0.2      59
#> 509         3   53.47697  interim         100       1, 1   0.2      60
#> 510         3   54.07651  interim         100       1, 1   0.2      61
#> 511         3   55.83878  interim         100       1, 1   0.2      62
#> 512         3   56.41014  interim         100       1, 1   0.2      63
#> 513         3   57.40339  interim         100       1, 1   0.2      64
#> 514         3   57.86333  interim         100       1, 1   0.2      65
#> 515         3   58.01604  interim         100       1, 1   0.2      66
#> 516         3   58.64184  interim         100       1, 1   0.2      67
#> 517         3   61.80970  interim         100       1, 1   0.2      68
#> 518         3   61.82024  interim         100       1, 1   0.2      69
#> 519         3   63.61490  interim         100       1, 1   0.2      70
#> 520         3   63.86728  interim         100       1, 1   0.2      70
#> 521         3   65.07797  interim         100       1, 1   0.2      71
#> 522         3   65.68095  interim         100       1, 1   0.2      72
#> 523         3   65.99680  interim         100       1, 1   0.2      73
#> 524         3   66.30091  interim         100       1, 1   0.2      74
#> 525         3   67.98914  interim         100       1, 1   0.2      75
#> 526         3   68.83942  interim         100       1, 1   0.2      76
#> 527         3   69.11459  interim         100       1, 1   0.2      77
#> 528         3   69.86589  interim         100       1, 1   0.2      78
#> 529         3   71.10411  interim         100       1, 1   0.2      79
#> 530         3   71.34250  interim         100       1, 1   0.2      80
#> 531         3   71.42586  interim         100       1, 1   0.2      81
#> 532         3   72.30656  interim         100       1, 1   0.2      82
#> 533         3   72.76776  interim         100       1, 1   0.2      83
#> 534         3   74.22738  interim         100       1, 1   0.2      84
#> 535         3   77.54871  interim         100       1, 1   0.2      85
#> 536         3   77.95592  interim         100       1, 1   0.2      86
#> 537         3   79.51993  interim         100       1, 1   0.2      87
#> 538         3   80.05824  interim         100       1, 1   0.2      88
#> 539         3   81.85533  interim         100       1, 1   0.2      89
#> 540         3   84.52151  interim         100       1, 1   0.2      90
#> 541         3   84.66379  interim         100       1, 1   0.2      91
#> 542         3   85.02198  interim         100       1, 1   0.2      92
#> 543         3   87.53058  interim         100       1, 1   0.2      93
#> 544         3   88.10133  interim         100       1, 1   0.2      94
#> 545         3   89.59937  interim         100       1, 1   0.2      95
#> 546         3   91.39456  interim         100       1, 1   0.2      96
#> 547         3   93.36668  interim         100       1, 1   0.2      97
#> 548         3   94.49027  interim         100       1, 1   0.2      98
#> 549         3   96.32723  interim         100       1, 1   0.2      99
#> 550         3   96.44952  interim         100       1, 1   0.2     100
#> 551         3   96.44952    final         100       1, 1   0.2     100
#> 552         3  641.30535  interim         100       1, 1   0.2     100
#> 553         3  641.30535    final         100       1, 1   0.2     100
#> 554         3  642.31945  interim         100       1, 1   0.2     100
#> 555         3  642.31945    final         100       1, 1   0.2     100
#> 556         3  674.93478  interim         100       1, 1   0.2     100
#> 557         3  674.93478    final         100       1, 1   0.2     100
#> 558         3  769.22693  interim         100       1, 1   0.2     100
#> 559         3  769.22693    final         100       1, 1   0.2     100
#> 560         3  842.08596  interim         100       1, 1   0.2     100
#> 561         3  842.08596    final         100       1, 1   0.2     100
#> 562         3  875.35735  interim         100       1, 1   0.2     100
#> 563         3  875.35735    final         100       1, 1   0.2     100
#> 564         3  907.68533  interim         100       1, 1   0.2     100
#> 565         3  907.68533    final         100       1, 1   0.2     100
#> 566         3  997.34185  interim         100       1, 1   0.2     100
#> 567         3  997.34185    final         100       1, 1   0.2     100
#> 568         3 1072.05718  interim         100       1, 1   0.2     100
#> 569         3 1072.05718    final         100       1, 1   0.2     100
#> 570         3 1163.39035  interim         100       1, 1   0.2     100
#> 571         3 1163.39035    final         100       1, 1   0.2     100
#> 572         3 1191.59587  interim         100       1, 1   0.2     100
#> 573         3 1191.59587    final         100       1, 1   0.2     100
#> 574         3 1562.07009  interim         100       1, 1   0.2     100
#> 575         3 1562.07009    final         100       1, 1   0.2     100
#> 576         3 1899.90611  interim         100       1, 1   0.2     100
#> 577         3 1899.90611    final         100       1, 1   0.2     100
#> 578         3 1940.35632  interim         100       1, 1   0.2     100
#> 579         3 1940.35632    final         100       1, 1   0.2     100
#> 580         3 2154.49290  interim         100       1, 1   0.2     100
#> 581         3 2154.49290    final         100       1, 1   0.2     100
#> 582         3 2447.37535  interim         100       1, 1   0.2     100
#> 583         3 2447.37535    final         100       1, 1   0.2     100
#> 584         3 2822.93260  interim         100       1, 1   0.2     100
#> 585         3 2822.93260    final         100       1, 1   0.2     100
#> 586         3 2941.35975  interim         100       1, 1   0.2     100
#> 587         3 2941.35975    final         100       1, 1   0.2     100
#> 588         3 3110.94827  interim         100       1, 1   0.2     100
#> 589         3 3110.94827    final         100       1, 1   0.2     100
#> 590         3 3148.98925  interim         100       1, 1   0.2     100
#> 591         3 3148.98925    final         100       1, 1   0.2     100
#> 592         3 3217.53725  interim         100       1, 1   0.2     100
#> 593         3 3217.53725    final         100       1, 1   0.2     100
#> 594         3 3337.18543  interim         100       1, 1   0.2     100
#> 595         3 3337.18543    final         100       1, 1   0.2     100
#> 596         3 3367.90075  interim         100       1, 1   0.2     100
#> 597         3 3367.90075    final         100       1, 1   0.2     100
#> 598         3 3529.54127  interim         100       1, 1   0.2     100
#> 599         3 3529.54127    final         100       1, 1   0.2     100
#> 600         3 3561.12345  interim         100       1, 1   0.2     100
#> 601         3 3561.12345    final         100       1, 1   0.2     100
#> 602         3 3709.40028  interim         100       1, 1   0.2     100
#> 603         3 3709.40028    final         100       1, 1   0.2     100
#> 604         3 3840.98646  interim         100       1, 1   0.2     100
#> 605         3 3840.98646    final         100       1, 1   0.2     100
#> 606         3 3981.35337  interim         100       1, 1   0.2     100
#> 607         3 3981.35337    final         100       1, 1   0.2     100
#> 608         3 4027.76788  interim         100       1, 1   0.2     100
#> 609         3 4027.76788    final         100       1, 1   0.2     100
#> 610         3 4160.72530  interim         100       1, 1   0.2     100
#> 611         3 4160.72530    final         100       1, 1   0.2     100
#> 612         3 4186.89686  interim         100       1, 1   0.2     100
#> 613         3 4186.89686    final         100       1, 1   0.2     100
#> 614         3 4399.55299  interim         100       1, 1   0.2     100
#> 615         3 4399.55299    final         100       1, 1   0.2     100
#> 616         3 4557.13615  interim         100       1, 1   0.2     100
#> 617         3 4557.13615    final         100       1, 1   0.2     100
#> 618         3 4638.33157  interim         100       1, 1   0.2     100
#> 619         3 4638.33157    final         100       1, 1   0.2     100
#> 620         3 4639.82467  interim         100       1, 1   0.2     100
#> 621         3 4639.82467    final         100       1, 1   0.2     100
#> 622         3 4671.16697  interim         100       1, 1   0.2     100
#> 623         3 4671.16697    final         100       1, 1   0.2     100
#> 624         3 4728.99685  interim         100       1, 1   0.2     100
#> 625         3 4728.99685    final         100       1, 1   0.2     100
#> 626         3 4770.76598  interim         100       1, 1   0.2     100
#> 627         3 4770.76598    final         100       1, 1   0.2     100
#> 628         3 4982.05892  interim         100       1, 1   0.2     100
#> 629         3 4982.05892    final         100       1, 1   0.2     100
#> 630         3 4984.32119  interim         100       1, 1   0.2     100
#> 631         3 4984.32119    final         100       1, 1   0.2     100
#> 632         3 5183.69330  interim         100       1, 1   0.2     100
#> 633         3 5183.69330    final         100       1, 1   0.2     100
#> 634         3 5246.51642  interim         100       1, 1   0.2     100
#> 635         3 5246.51642    final         100       1, 1   0.2     100
#> 636         3 5277.27244  interim         100       1, 1   0.2     100
#> 637         3 5277.27244    final         100       1, 1   0.2     100
#> 638         3 5387.88795  interim         100       1, 1   0.2     100
#> 639         3 5387.88795    final         100       1, 1   0.2     100
#> 640         3 5407.20693  interim         100       1, 1   0.2     100
#> 641         3 5407.20693    final         100       1, 1   0.2     100
#> 642         3 5452.74780  interim         100       1, 1   0.2     100
#> 643         3 5452.74780    final         100       1, 1   0.2     100
#> 644         3 5468.71062  interim         100       1, 1   0.2     100
#> 645         3 5468.71062    final         100       1, 1   0.2     100
#> 646         3 5542.58750  interim         100       1, 1   0.2     100
#> 647         3 5542.58750    final         100       1, 1   0.2     100
#> 648         3 5763.47315  interim         100       1, 1   0.2     100
#> 649         3 5763.47315    final         100       1, 1   0.2     100
#> 650         3 5771.96148  interim         100       1, 1   0.2     100
#> 651         3 5771.96148    final         100       1, 1   0.2     100
#> 652         3 5984.73519  interim         100       1, 1   0.2     100
#> 653         3 5984.73519    final         100       1, 1   0.2     100
#> 654         3 6004.08873  interim         100       1, 1   0.2     100
#> 655         3 6004.08873    final         100       1, 1   0.2     100
#> 656         3 6300.91345  interim         100       1, 1   0.2     100
#> 657         3 6300.91345    final         100       1, 1   0.2     100
#> 658         3 6311.73592  interim         100       1, 1   0.2     100
#> 659         3 6311.73592    final         100       1, 1   0.2     100
#> 660         3 6313.31644  interim         100       1, 1   0.2     100
#> 661         3 6313.31644    final         100       1, 1   0.2     100
#> 662         3 6417.45507  interim         100       1, 1   0.2     100
#> 663         3 6417.45507    final         100       1, 1   0.2     100
#> 664         3 6569.04139  interim         100       1, 1   0.2     100
#> 665         3 6569.04139    final         100       1, 1   0.2     100
#> 666         3 6679.38170  interim         100       1, 1   0.2     100
#> 667         3 6679.38170    final         100       1, 1   0.2     100
#> 668         3 6809.27436  interim         100       1, 1   0.2     100
#> 669         3 6809.27436    final         100       1, 1   0.2     100
#> 670         3 6813.54794  interim         100       1, 1   0.2     100
#> 671         3 6813.54794    final         100       1, 1   0.2     100
#> 672         3 6849.44519  interim         100       1, 1   0.2     100
#> 673         3 6849.44519    final         100       1, 1   0.2     100
#> 674         3 6970.02494  interim         100       1, 1   0.2     100
#> 675         3 6970.02494    final         100       1, 1   0.2     100
#> 676         3 6987.69901  interim         100       1, 1   0.2     100
#> 677         3 6987.69901    final         100       1, 1   0.2     100
#> 678         3 7057.59733  interim         100       1, 1   0.2     100
#> 679         3 7057.59733    final         100       1, 1   0.2     100
#> 680         3 7120.24351  interim         100       1, 1   0.2     100
#> 681         3 7120.24351    final         100       1, 1   0.2     100
#> 682         3 7148.95843  interim         100       1, 1   0.2     100
#> 683         3 7148.95843    final         100       1, 1   0.2     100
#> 684         3 7367.40178  interim         100       1, 1   0.2     100
#> 685         3 7367.40178    final         100       1, 1   0.2     100
#> 686         3 7413.14504  interim         100       1, 1   0.2     100
#> 687         3 7413.14504    final         100       1, 1   0.2     100
#> 688         3 7447.56728  interim         100       1, 1   0.2     100
#> 689         3 7447.56728    final         100       1, 1   0.2     100
#> 690         3 7511.52748  interim         100       1, 1   0.2     100
#> 691         3 7511.52748    final         100       1, 1   0.2     100
#> 692         3 7545.62878  interim         100       1, 1   0.2     100
#> 693         3 7545.62878    final         100       1, 1   0.2     100
#> 694         3 7575.01451  interim         100       1, 1   0.2     100
#> 695         3 7575.01451    final         100       1, 1   0.2     100
#> 696         3 7653.02453  interim         100       1, 1   0.2     100
#> 697         3 7653.02453    final         100       1, 1   0.2     100
#> 698         3 7718.38987  interim         100       1, 1   0.2     100
#> 699         3 7718.38987    final         100       1, 1   0.2     100
#> 700         3 7742.26552  interim         100       1, 1   0.2     100
#> 701         3 7742.26552    final         100       1, 1   0.2     100
#> 702         3 7810.47626  interim         100       1, 1   0.2     100
#> 703         3 7810.47626    final         100       1, 1   0.2     100
#> 704         3 7929.15282  interim         100       1, 1   0.2     100
#> 705         3 7929.15282    final         100       1, 1   0.2     100
#> 706         3 8138.48451  interim         100       1, 1   0.2     100
#> 707         3 8138.48451    final         100       1, 1   0.2     100
#> 708         3 8213.91817  interim         100       1, 1   0.2     100
#> 709         3 8213.91817    final         100       1, 1   0.2     100
#> 710         3 8356.04399  interim         100       1, 1   0.2     100
#> 711         3 8356.04399    final         100       1, 1   0.2     100
#> 712         3 8394.89321  interim         100       1, 1   0.2     100
#> 713         3 8394.89321    final         100       1, 1   0.2     100
#> 714         3 8567.68352  interim         100       1, 1   0.2     100
#> 715         3 8567.68352    final         100       1, 1   0.2     100
#> 716         3 8612.32526  interim         100       1, 1   0.2     100
#> 717         3 8612.32526    final         100       1, 1   0.2     100
#> 718         3 8683.38045  interim         100       1, 1   0.2     100
#> 719         3 8683.38045    final         100       1, 1   0.2     100
#> 720         3 8693.11929  interim         100       1, 1   0.2     100
#> 721         3 8693.11929    final         100       1, 1   0.2     100
#> 722         3 8832.87024  interim         100       1, 1   0.2     100
#> 723         3 8832.87024    final         100       1, 1   0.2     100
#> 724         3 8900.00150  interim         100       1, 1   0.2     100
#> 725         3 8900.00150    final         100       1, 1   0.2     100
#> 726         3 8988.58802  interim         100       1, 1   0.2     100
#> 727         3 8988.58802    final         100       1, 1   0.2     100
#> 728         3 9110.43423  interim         100       1, 1   0.2     100
#> 729         3 9110.43423    final         100       1, 1   0.2     100
#> 730         3 9272.97598  interim         100       1, 1   0.2     100
#> 731         3 9272.97598    final         100       1, 1   0.2     100
#> 732         3 9283.41623  interim         100       1, 1   0.2     100
#> 733         3 9283.41623    final         100       1, 1   0.2     100
#> 734         3 9313.77786  interim         100       1, 1   0.2     100
#> 735         3 9313.77786    final         100       1, 1   0.2     100
#> 736         3 9394.88550  interim         100       1, 1   0.2     100
#> 737         3 9394.88550    final         100       1, 1   0.2     100
#> 738         3 9562.00349  interim         100       1, 1   0.2     100
#> 739         3 9562.00349    final         100       1, 1   0.2     100
#> 740         3 9594.28157  interim         100       1, 1   0.2     100
#> 741         3 9594.28157    final         100       1, 1   0.2     100
#> 742         3 9701.40004  interim         100       1, 1   0.2     100
#> 743         3 9701.40004    final         100       1, 1   0.2     100
#> 744         3 9706.49654  interim         100       1, 1   0.2     100
#> 745         3 9706.49654    final         100       1, 1   0.2     100
#> 746         3 9730.28383  interim         100       1, 1   0.2     100
#> 747         3 9730.28383    final         100       1, 1   0.2     100
#>        mean_pbo    mean_trt    p_value
#> 1   -0.15458732  0.42624025 0.01626231
#> 2   -0.15458732  0.37749229 0.02687530
#> 3   -0.15458732  0.37414045 0.02464318
#> 4   -0.16888919  0.37414045 0.01925765
#> 5   -0.16888919  0.44412744 0.01055039
#> 6   -0.16888919  0.41689247 0.01299920
#> 7   -0.13071314  0.41689247 0.01951807
#> 8   -0.11762134  0.41689247 0.02072547
#> 9   -0.11762134  0.42515236 0.01671901
#> 10  -0.07328431  0.42515236 0.02812241
#> 11  -0.07328431  0.42586599 0.02511221
#> 12  -0.09211570  0.42586599 0.01894600
#> 13  -0.05942002  0.42586599 0.02703900
#> 14  -0.05942002  0.36181621 0.06006327
#> 15  -0.01578277  0.36181621 0.09297339
#> 16  -0.01578277  0.31269790 0.14545451
#> 17  -0.01578277  0.31000738 0.14168045
#> 18  -0.03113731  0.31000738 0.12017204
#> 19  -0.04401287  0.31000738 0.10338143
#> 20  -0.04401287  0.34102245 0.07506579
#> 21  -0.02048135  0.34102245 0.09231136
#> 22  -0.01156076  0.34102245 0.09695881
#> 23  -0.04666896  0.34102245 0.06943028
#> 24  -0.04666896  0.36277489 0.05270111
#> 25  -0.05559799  0.36277489 0.04574890
#> 26  -0.09112265  0.36277489 0.03134163
#> 27  -0.09112265  0.32562866 0.04742676
#> 28  -0.08152389  0.32562866 0.05054836
#> 29  -0.08152389  0.25838981 0.11383775
#> 30  -0.10004798  0.25838981 0.09410191
#> 31  -0.10004798  0.24784925 0.09861576
#> 32  -0.10004798  0.30361838 0.06094877
#> 33  -0.09191048  0.30361838 0.06417478
#> 34  -0.09191048  0.30842468 0.05697569
#> 35  -0.09191048  0.30510754 0.05513630
#> 36  -0.05896468  0.30510754 0.07931298
#> 37  -0.05896468  0.30399454 0.07578260
#> 38  -0.06163035  0.30399454 0.07124254
#> 39  -0.06163035  0.33522514 0.05011246
#> 40  -0.06163035  0.35355699 0.03867769
#> 41  -0.06163035  0.35867032 0.03397010
#> 42  -0.06163035  0.38682777 0.02354243
#> 43  -0.05545868  0.38682777 0.02428096
#> 44  -0.01864018  0.38682777 0.04026614
#> 45  -0.01864018  0.36436087 0.05098385
#> 46  -0.02804583  0.36436087 0.04408675
#> 47  -0.02804583  0.36011971 0.04373186
#> 48  -0.02711106  0.36011971 0.04247915
#> 49  -0.02711106  0.35801094 0.04107039
#> 50  -0.02475451  0.35801094 0.04064238
#> 51  -0.02475451  0.32833062 0.05877739
#> 52  -0.02475451  0.32833062 0.05877739
#> 53  -0.02475451  0.32833062 0.05877739
#> 54  -0.02475451  0.32833062 0.05877739
#> 55  -0.02475451  0.32833062 0.05877739
#> 56  -0.02475451  0.32833062 0.05877739
#> 57  -0.02475451  0.32833062 0.05877739
#> 58  -0.02475451  0.32833062 0.05877739
#> 59  -0.02475451  0.32833062 0.05877739
#> 60  -0.02475451  0.32833062 0.05877739
#> 61  -0.02475451  0.32833062 0.05877739
#> 62  -0.02475451  0.32833062 0.05877739
#> 63  -0.02475451  0.32833062 0.05877739
#> 64  -0.02475451  0.32833062 0.05877739
#> 65  -0.02475451  0.32833062 0.05877739
#> 66  -0.02475451  0.32833062 0.05877739
#> 67  -0.02475451  0.32833062 0.05877739
#> 68  -0.02475451  0.32833062 0.05877739
#> 69  -0.02475451  0.32833062 0.05877739
#> 70  -0.02475451  0.32833062 0.05877739
#> 71  -0.02475451  0.32833062 0.05877739
#> 72  -0.02475451  0.32833062 0.05877739
#> 73  -0.02475451  0.32833062 0.05877739
#> 74  -0.02475451  0.32833062 0.05877739
#> 75  -0.02475451  0.32833062 0.05877739
#> 76  -0.02475451  0.32833062 0.05877739
#> 77  -0.02475451  0.32833062 0.05877739
#> 78  -0.02475451  0.32833062 0.05877739
#> 79  -0.02475451  0.32833062 0.05877739
#> 80  -0.02475451  0.32833062 0.05877739
#> 81  -0.02475451  0.32833062 0.05877739
#> 82  -0.02475451  0.32833062 0.05877739
#> 83  -0.02475451  0.32833062 0.05877739
#> 84  -0.02475451  0.32833062 0.05877739
#> 85  -0.02475451  0.32833062 0.05877739
#> 86  -0.02475451  0.32833062 0.05877739
#> 87  -0.02475451  0.32833062 0.05877739
#> 88  -0.02475451  0.32833062 0.05877739
#> 89  -0.02475451  0.32833062 0.05877739
#> 90  -0.02475451  0.32833062 0.05877739
#> 91  -0.02475451  0.32833062 0.05877739
#> 92  -0.02475451  0.32833062 0.05877739
#> 93  -0.02475451  0.32833062 0.05877739
#> 94  -0.02475451  0.32833062 0.05877739
#> 95  -0.02475451  0.32833062 0.05877739
#> 96  -0.02475451  0.32833062 0.05877739
#> 97  -0.02475451  0.32833062 0.05877739
#> 98  -0.02475451  0.32833062 0.05877739
#> 99  -0.02475451  0.32833062 0.05877739
#> 100 -0.02475451  0.32833062 0.05877739
#> 101 -0.02475451  0.32833062 0.05877739
#> 102 -0.02475451  0.32833062 0.05877739
#> 103 -0.02475451  0.32833062 0.05877739
#> 104 -0.02475451  0.32833062 0.05877739
#> 105 -0.02475451  0.32833062 0.05877739
#> 106 -0.02475451  0.32833062 0.05877739
#> 107 -0.02475451  0.32833062 0.05877739
#> 108 -0.02475451  0.32833062 0.05877739
#> 109 -0.02475451  0.32833062 0.05877739
#> 110 -0.02475451  0.32833062 0.05877739
#> 111 -0.02475451  0.32833062 0.05877739
#> 112 -0.02475451  0.32833062 0.05877739
#> 113 -0.02475451  0.32833062 0.05877739
#> 114 -0.02475451  0.32833062 0.05877739
#> 115 -0.02475451  0.32833062 0.05877739
#> 116 -0.02475451  0.32833062 0.05877739
#> 117 -0.02475451  0.32833062 0.05877739
#> 118 -0.02475451  0.32833062 0.05877739
#> 119 -0.02475451  0.32833062 0.05877739
#> 120 -0.02475451  0.32833062 0.05877739
#> 121 -0.02475451  0.32833062 0.05877739
#> 122 -0.02475451  0.32833062 0.05877739
#> 123 -0.02475451  0.32833062 0.05877739
#> 124 -0.02475451  0.32833062 0.05877739
#> 125 -0.02475451  0.32833062 0.05877739
#> 126 -0.02475451  0.32833062 0.05877739
#> 127 -0.02475451  0.32833062 0.05877739
#> 128 -0.02475451  0.32833062 0.05877739
#> 129 -0.02475451  0.32833062 0.05877739
#> 130 -0.02475451  0.32833062 0.05877739
#> 131 -0.02475451  0.32833062 0.05877739
#> 132 -0.02475451  0.32833062 0.05877739
#> 133 -0.02475451  0.32833062 0.05877739
#> 134 -0.02475451  0.32833062 0.05877739
#> 135 -0.02475451  0.32833062 0.05877739
#> 136 -0.02475451  0.32833062 0.05877739
#> 137 -0.02475451  0.32833062 0.05877739
#> 138 -0.02475451  0.32833062 0.05877739
#> 139 -0.02475451  0.32833062 0.05877739
#> 140 -0.02475451  0.32833062 0.05877739
#> 141 -0.02475451  0.32833062 0.05877739
#> 142 -0.02475451  0.32833062 0.05877739
#> 143 -0.02475451  0.32833062 0.05877739
#> 144 -0.02475451  0.32833062 0.05877739
#> 145 -0.02475451  0.32833062 0.05877739
#> 146 -0.02475451  0.32833062 0.05877739
#> 147 -0.02475451  0.32833062 0.05877739
#> 148 -0.02475451  0.32833062 0.05877739
#> 149 -0.02475451  0.32833062 0.05877739
#> 150 -0.02475451  0.32833062 0.05877739
#> 151 -0.02475451  0.32833062 0.05877739
#> 152 -0.02475451  0.32833062 0.05877739
#> 153 -0.02475451  0.32833062 0.05877739
#> 154 -0.02475451  0.32833062 0.05877739
#> 155 -0.02475451  0.32833062 0.05877739
#> 156 -0.02475451  0.32833062 0.05877739
#> 157 -0.02475451  0.32833062 0.05877739
#> 158 -0.02475451  0.32833062 0.05877739
#> 159 -0.02475451  0.32833062 0.05877739
#> 160 -0.02475451  0.32833062 0.05877739
#> 161 -0.02475451  0.32833062 0.05877739
#> 162 -0.02475451  0.32833062 0.05877739
#> 163 -0.02475451  0.32833062 0.05877739
#> 164 -0.02475451  0.32833062 0.05877739
#> 165 -0.02475451  0.32833062 0.05877739
#> 166 -0.02475451  0.32833062 0.05877739
#> 167 -0.02475451  0.32833062 0.05877739
#> 168 -0.02475451  0.32833062 0.05877739
#> 169 -0.02475451  0.32833062 0.05877739
#> 170 -0.02475451  0.32833062 0.05877739
#> 171 -0.02475451  0.32833062 0.05877739
#> 172 -0.02475451  0.32833062 0.05877739
#> 173 -0.02475451  0.32833062 0.05877739
#> 174 -0.02475451  0.32833062 0.05877739
#> 175 -0.02475451  0.32833062 0.05877739
#> 176 -0.02475451  0.32833062 0.05877739
#> 177 -0.02475451  0.32833062 0.05877739
#> 178 -0.02475451  0.32833062 0.05877739
#> 179 -0.02475451  0.32833062 0.05877739
#> 180 -0.02475451  0.32833062 0.05877739
#> 181 -0.02475451  0.32833062 0.05877739
#> 182 -0.02475451  0.32833062 0.05877739
#> 183 -0.02475451  0.32833062 0.05877739
#> 184 -0.02475451  0.32833062 0.05877739
#> 185 -0.02475451  0.32833062 0.05877739
#> 186 -0.02475451  0.32833062 0.05877739
#> 187 -0.02475451  0.32833062 0.05877739
#> 188 -0.02475451  0.32833062 0.05877739
#> 189 -0.02475451  0.32833062 0.05877739
#> 190 -0.02475451  0.32833062 0.05877739
#> 191 -0.02475451  0.32833062 0.05877739
#> 192 -0.02475451  0.32833062 0.05877739
#> 193 -0.02475451  0.32833062 0.05877739
#> 194 -0.02475451  0.32833062 0.05877739
#> 195 -0.02475451  0.32833062 0.05877739
#> 196 -0.02475451  0.32833062 0.05877739
#> 197 -0.02475451  0.32833062 0.05877739
#> 198 -0.02475451  0.32833062 0.05877739
#> 199 -0.02475451  0.32833062 0.05877739
#> 200 -0.02475451  0.32833062 0.05877739
#> 201 -0.02475451  0.32833062 0.05877739
#> 202 -0.02475451  0.32833062 0.05877739
#> 203 -0.02475451  0.32833062 0.05877739
#> 204 -0.02475451  0.32833062 0.05877739
#> 205 -0.02475451  0.32833062 0.05877739
#> 206 -0.02475451  0.32833062 0.05877739
#> 207 -0.02475451  0.32833062 0.05877739
#> 208 -0.02475451  0.32833062 0.05877739
#> 209 -0.02475451  0.32833062 0.05877739
#> 210 -0.02475451  0.32833062 0.05877739
#> 211 -0.02475451  0.32833062 0.05877739
#> 212 -0.02475451  0.32833062 0.05877739
#> 213 -0.02475451  0.32833062 0.05877739
#> 214 -0.02475451  0.32833062 0.05877739
#> 215 -0.02475451  0.32833062 0.05877739
#> 216 -0.02475451  0.32833062 0.05877739
#> 217 -0.02475451  0.32833062 0.05877739
#> 218 -0.02475451  0.32833062 0.05877739
#> 219 -0.02475451  0.32833062 0.05877739
#> 220 -0.02475451  0.32833062 0.05877739
#> 221 -0.02475451  0.32833062 0.05877739
#> 222 -0.02475451  0.32833062 0.05877739
#> 223 -0.02475451  0.32833062 0.05877739
#> 224 -0.02475451  0.32833062 0.05877739
#> 225 -0.02475451  0.32833062 0.05877739
#> 226 -0.02475451  0.32833062 0.05877739
#> 227 -0.02475451  0.32833062 0.05877739
#> 228 -0.02475451  0.32833062 0.05877739
#> 229 -0.02475451  0.32833062 0.05877739
#> 230 -0.02475451  0.32833062 0.05877739
#> 231 -0.02475451  0.32833062 0.05877739
#> 232 -0.02475451  0.32833062 0.05877739
#> 233 -0.02475451  0.32833062 0.05877739
#> 234 -0.02475451  0.32833062 0.05877739
#> 235 -0.02475451  0.32833062 0.05877739
#> 236 -0.02475451  0.32833062 0.05877739
#> 237 -0.02475451  0.32833062 0.05877739
#> 238 -0.02475451  0.32833062 0.05877739
#> 239 -0.02475451  0.32833062 0.05877739
#> 240 -0.02475451  0.32833062 0.05877739
#> 241 -0.02475451  0.32833062 0.05877739
#> 242 -0.02475451  0.32833062 0.05877739
#> 243 -0.02475451  0.32833062 0.05877739
#> 244 -0.02475451  0.32833062 0.05877739
#> 245 -0.02475451  0.32833062 0.05877739
#> 246 -0.02475451  0.32833062 0.05877739
#> 247 -0.04934773  0.02698932 0.80070480
#> 248 -0.04934773  0.09869826 0.62851043
#> 249 -0.04934773  0.11210911 0.59179596
#> 250 -0.04934773  0.11288434 0.58400256
#> 251 -0.02989248  0.11288434 0.62262264
#> 252 -0.06299402  0.11288434 0.53875938
#> 253 -0.06299402  0.07251811 0.63354879
#> 254 -0.01933395  0.07251811 0.74450040
#> 255 -0.01933395  0.11574167 0.63073193
#> 256 -0.02635203  0.11574167 0.60633935
#> 257 -0.02635203  0.11079665 0.61345738
#> 258 -0.02635203  0.13454551 0.54937007
#> 259 -0.07676065  0.13454551 0.43183720
#> 260 -0.07676065  0.19605435 0.31704177
#> 261 -0.05182129  0.19605435 0.35650113
#> 262 -0.05182129  0.20246975 0.33758030
#> 263 -0.03011311  0.20246975 0.37357764
#> 264 -0.03011311  0.22276947 0.32845013
#> 265 -0.03011311  0.23749916 0.29578768
#> 266 -0.08780447  0.23749916 0.20867554
#> 267 -0.04345569  0.23749916 0.27605426
#> 268 -0.03912615  0.23749916 0.27572826
#> 269 -0.04133790  0.23749916 0.26464655
#> 270 -0.09933677  0.23749916 0.18424889
#> 271 -0.09933677  0.26489751 0.14878323
#> 272 -0.11566719  0.26489751 0.12692598
#> 273 -0.11566719  0.28739243 0.10345860
#> 274 -0.12875964  0.28739243 0.08884740
#> 275 -0.11680052  0.28739243 0.09420134
#> 276 -0.11680052  0.27597666 0.09981568
#> 277 -0.11680052  0.26254818 0.10809972
#> 278 -0.08261549  0.26254818 0.14215610
#> 279 -0.04608380  0.26254818 0.18850760
#> 280 -0.04608380  0.25045505 0.20167221
#> 281 -0.03159353  0.25045505 0.21934224
#> 282 -0.03159353  0.18505220 0.35849950
#> 283 -0.03159353  0.18790621 0.34666303
#> 284 -0.03943525  0.18790621 0.32432592
#> 285 -0.03943525  0.18485079 0.32537841
#> 286 -0.04013892  0.18485079 0.31824878
#> 287 -0.06017972  0.18485079 0.27392780
#> 288 -0.06619135  0.18485079 0.25750237
#> 289 -0.06619135  0.17962848 0.26208778
#> 290 -0.06619135  0.15292571 0.31545548
#> 291 -0.06619135  0.12543989 0.37840896
#> 292 -0.07363654  0.12543989 0.35544048
#> 293 -0.07363654  0.11737571 0.37041223
#> 294 -0.09430768  0.11737571 0.31843280
#> 295 -0.09430768  0.13052493 0.28540743
#> 296 -0.08461093  0.13052493 0.30217147
#> 297 -0.05831826  0.13052493 0.36397787
#> 298 -0.05831826  0.13052493 0.36397787
#> 299 -0.05831826  0.13052493 0.36397787
#> 300 -0.05831826  0.13052493 0.36397787
#> 301 -0.05831826  0.13052493 0.36397787
#> 302 -0.05831826  0.13052493 0.36397787
#> 303 -0.05831826  0.13052493 0.36397787
#> 304 -0.05831826  0.13052493 0.36397787
#> 305 -0.05831826  0.13052493 0.36397787
#> 306 -0.05831826  0.13052493 0.36397787
#> 307 -0.05831826  0.13052493 0.36397787
#> 308 -0.05831826  0.13052493 0.36397787
#> 309 -0.05831826  0.13052493 0.36397787
#> 310 -0.05831826  0.13052493 0.36397787
#> 311 -0.05831826  0.13052493 0.36397787
#> 312 -0.05831826  0.13052493 0.36397787
#> 313 -0.05831826  0.13052493 0.36397787
#> 314 -0.05831826  0.13052493 0.36397787
#> 315 -0.05831826  0.13052493 0.36397787
#> 316 -0.05831826  0.13052493 0.36397787
#> 317 -0.05831826  0.13052493 0.36397787
#> 318 -0.05831826  0.13052493 0.36397787
#> 319 -0.05831826  0.13052493 0.36397787
#> 320 -0.05831826  0.13052493 0.36397787
#> 321 -0.05831826  0.13052493 0.36397787
#> 322 -0.05831826  0.13052493 0.36397787
#> 323 -0.05831826  0.13052493 0.36397787
#> 324 -0.05831826  0.13052493 0.36397787
#> 325 -0.05831826  0.13052493 0.36397787
#> 326 -0.05831826  0.13052493 0.36397787
#> 327 -0.05831826  0.13052493 0.36397787
#> 328 -0.05831826  0.13052493 0.36397787
#> 329 -0.05831826  0.13052493 0.36397787
#> 330 -0.05831826  0.13052493 0.36397787
#> 331 -0.05831826  0.13052493 0.36397787
#> 332 -0.05831826  0.13052493 0.36397787
#> 333 -0.05831826  0.13052493 0.36397787
#> 334 -0.05831826  0.13052493 0.36397787
#> 335 -0.05831826  0.13052493 0.36397787
#> 336 -0.05831826  0.13052493 0.36397787
#> 337 -0.05831826  0.13052493 0.36397787
#> 338 -0.05831826  0.13052493 0.36397787
#> 339 -0.05831826  0.13052493 0.36397787
#> 340 -0.05831826  0.13052493 0.36397787
#> 341 -0.05831826  0.13052493 0.36397787
#> 342 -0.05831826  0.13052493 0.36397787
#> 343 -0.05831826  0.13052493 0.36397787
#> 344 -0.05831826  0.13052493 0.36397787
#> 345 -0.05831826  0.13052493 0.36397787
#> 346 -0.05831826  0.13052493 0.36397787
#> 347 -0.05831826  0.13052493 0.36397787
#> 348 -0.05831826  0.13052493 0.36397787
#> 349 -0.05831826  0.13052493 0.36397787
#> 350 -0.05831826  0.13052493 0.36397787
#> 351 -0.05831826  0.13052493 0.36397787
#> 352 -0.05831826  0.13052493 0.36397787
#> 353 -0.05831826  0.13052493 0.36397787
#> 354 -0.05831826  0.13052493 0.36397787
#> 355 -0.05831826  0.13052493 0.36397787
#> 356 -0.05831826  0.13052493 0.36397787
#> 357 -0.05831826  0.13052493 0.36397787
#> 358 -0.05831826  0.13052493 0.36397787
#> 359 -0.05831826  0.13052493 0.36397787
#> 360 -0.05831826  0.13052493 0.36397787
#> 361 -0.05831826  0.13052493 0.36397787
#> 362 -0.05831826  0.13052493 0.36397787
#> 363 -0.05831826  0.13052493 0.36397787
#> 364 -0.05831826  0.13052493 0.36397787
#> 365 -0.05831826  0.13052493 0.36397787
#> 366 -0.05831826  0.13052493 0.36397787
#> 367 -0.05831826  0.13052493 0.36397787
#> 368 -0.05831826  0.13052493 0.36397787
#> 369 -0.05831826  0.13052493 0.36397787
#> 370 -0.05831826  0.13052493 0.36397787
#> 371 -0.05831826  0.13052493 0.36397787
#> 372 -0.05831826  0.13052493 0.36397787
#> 373 -0.05831826  0.13052493 0.36397787
#> 374 -0.05831826  0.13052493 0.36397787
#> 375 -0.05831826  0.13052493 0.36397787
#> 376 -0.05831826  0.13052493 0.36397787
#> 377 -0.05831826  0.13052493 0.36397787
#> 378 -0.05831826  0.13052493 0.36397787
#> 379 -0.05831826  0.13052493 0.36397787
#> 380 -0.05831826  0.13052493 0.36397787
#> 381 -0.05831826  0.13052493 0.36397787
#> 382 -0.05831826  0.13052493 0.36397787
#> 383 -0.05831826  0.13052493 0.36397787
#> 384 -0.05831826  0.13052493 0.36397787
#> 385 -0.05831826  0.13052493 0.36397787
#> 386 -0.05831826  0.13052493 0.36397787
#> 387 -0.05831826  0.13052493 0.36397787
#> 388 -0.05831826  0.13052493 0.36397787
#> 389 -0.05831826  0.13052493 0.36397787
#> 390 -0.05831826  0.13052493 0.36397787
#> 391 -0.05831826  0.13052493 0.36397787
#> 392 -0.05831826  0.13052493 0.36397787
#> 393 -0.05831826  0.13052493 0.36397787
#> 394 -0.05831826  0.13052493 0.36397787
#> 395 -0.05831826  0.13052493 0.36397787
#> 396 -0.05831826  0.13052493 0.36397787
#> 397 -0.05831826  0.13052493 0.36397787
#> 398 -0.05831826  0.13052493 0.36397787
#> 399 -0.05831826  0.13052493 0.36397787
#> 400 -0.05831826  0.13052493 0.36397787
#> 401 -0.05831826  0.13052493 0.36397787
#> 402 -0.05831826  0.13052493 0.36397787
#> 403 -0.05831826  0.13052493 0.36397787
#> 404 -0.05831826  0.13052493 0.36397787
#> 405 -0.05831826  0.13052493 0.36397787
#> 406 -0.05831826  0.13052493 0.36397787
#> 407 -0.05831826  0.13052493 0.36397787
#> 408 -0.05831826  0.13052493 0.36397787
#> 409 -0.05831826  0.13052493 0.36397787
#> 410 -0.05831826  0.13052493 0.36397787
#> 411 -0.05831826  0.13052493 0.36397787
#> 412 -0.05831826  0.13052493 0.36397787
#> 413 -0.05831826  0.13052493 0.36397787
#> 414 -0.05831826  0.13052493 0.36397787
#> 415 -0.05831826  0.13052493 0.36397787
#> 416 -0.05831826  0.13052493 0.36397787
#> 417 -0.05831826  0.13052493 0.36397787
#> 418 -0.05831826  0.13052493 0.36397787
#> 419 -0.05831826  0.13052493 0.36397787
#> 420 -0.05831826  0.13052493 0.36397787
#> 421 -0.05831826  0.13052493 0.36397787
#> 422 -0.05831826  0.13052493 0.36397787
#> 423 -0.05831826  0.13052493 0.36397787
#> 424 -0.05831826  0.13052493 0.36397787
#> 425 -0.05831826  0.13052493 0.36397787
#> 426 -0.05831826  0.13052493 0.36397787
#> 427 -0.05831826  0.13052493 0.36397787
#> 428 -0.05831826  0.13052493 0.36397787
#> 429 -0.05831826  0.13052493 0.36397787
#> 430 -0.05831826  0.13052493 0.36397787
#> 431 -0.05831826  0.13052493 0.36397787
#> 432 -0.05831826  0.13052493 0.36397787
#> 433 -0.05831826  0.13052493 0.36397787
#> 434 -0.05831826  0.13052493 0.36397787
#> 435 -0.05831826  0.13052493 0.36397787
#> 436 -0.05831826  0.13052493 0.36397787
#> 437 -0.05831826  0.13052493 0.36397787
#> 438 -0.05831826  0.13052493 0.36397787
#> 439 -0.05831826  0.13052493 0.36397787
#> 440 -0.05831826  0.13052493 0.36397787
#> 441 -0.05831826  0.13052493 0.36397787
#> 442 -0.05831826  0.13052493 0.36397787
#> 443 -0.05831826  0.13052493 0.36397787
#> 444 -0.05831826  0.13052493 0.36397787
#> 445 -0.05831826  0.13052493 0.36397787
#> 446 -0.05831826  0.13052493 0.36397787
#> 447 -0.05831826  0.13052493 0.36397787
#> 448 -0.05831826  0.13052493 0.36397787
#> 449 -0.05831826  0.13052493 0.36397787
#> 450 -0.05831826  0.13052493 0.36397787
#> 451 -0.05831826  0.13052493 0.36397787
#> 452 -0.05831826  0.13052493 0.36397787
#> 453 -0.05831826  0.13052493 0.36397787
#> 454 -0.05831826  0.13052493 0.36397787
#> 455 -0.05831826  0.13052493 0.36397787
#> 456 -0.05831826  0.13052493 0.36397787
#> 457 -0.05831826  0.13052493 0.36397787
#> 458 -0.05831826  0.13052493 0.36397787
#> 459 -0.05831826  0.13052493 0.36397787
#> 460 -0.05831826  0.13052493 0.36397787
#> 461 -0.05831826  0.13052493 0.36397787
#> 462 -0.05831826  0.13052493 0.36397787
#> 463 -0.05831826  0.13052493 0.36397787
#> 464 -0.05831826  0.13052493 0.36397787
#> 465 -0.05831826  0.13052493 0.36397787
#> 466 -0.05831826  0.13052493 0.36397787
#> 467 -0.05831826  0.13052493 0.36397787
#> 468 -0.05831826  0.13052493 0.36397787
#> 469 -0.05831826  0.13052493 0.36397787
#> 470 -0.05831826  0.13052493 0.36397787
#> 471 -0.05831826  0.13052493 0.36397787
#> 472 -0.05831826  0.13052493 0.36397787
#> 473 -0.05831826  0.13052493 0.36397787
#> 474 -0.05831826  0.13052493 0.36397787
#> 475 -0.05831826  0.13052493 0.36397787
#> 476 -0.05831826  0.13052493 0.36397787
#> 477 -0.05831826  0.13052493 0.36397787
#> 478 -0.05831826  0.13052493 0.36397787
#> 479 -0.05831826  0.13052493 0.36397787
#> 480 -0.05831826  0.13052493 0.36397787
#> 481 -0.05831826  0.13052493 0.36397787
#> 482 -0.05831826  0.13052493 0.36397787
#> 483 -0.05831826  0.13052493 0.36397787
#> 484 -0.05831826  0.13052493 0.36397787
#> 485 -0.05831826  0.13052493 0.36397787
#> 486 -0.05831826  0.13052493 0.36397787
#> 487 -0.05831826  0.13052493 0.36397787
#> 488 -0.05831826  0.13052493 0.36397787
#> 489 -0.05831826  0.13052493 0.36397787
#> 490 -0.05831826  0.13052493 0.36397787
#> 491 -0.05831826  0.13052493 0.36397787
#> 492 -0.05831826  0.13052493 0.36397787
#> 493 -0.05831826  0.13052493 0.36397787
#> 494 -0.05831826  0.13052493 0.36397787
#> 495 -0.05831826  0.13052493 0.36397787
#> 496 -0.05831826  0.13052493 0.36397787
#> 497 -0.05831826  0.13052493 0.36397787
#> 498 -0.05831826  0.13052493 0.36397787
#> 499 -0.14213989 -0.13708169 0.98578393
#> 500 -0.14213989 -0.07810005 0.82112099
#> 501 -0.14213989 -0.02455509 0.67661747
#> 502 -0.11274408 -0.02455509 0.75198759
#> 503 -0.12992344 -0.02455509 0.70228463
#> 504 -0.15160166 -0.02455509 0.64160245
#> 505 -0.13315922 -0.02455509 0.68761210
#> 506 -0.18229742 -0.02455509 0.56149801
#> 507 -0.23831953 -0.02455509 0.43650499
#> 508 -0.23831953 -0.02411999 0.42441575
#> 509 -0.23831953 -0.02261879 0.41071591
#> 510 -0.22585245 -0.02261879 0.43282312
#> 511 -0.27375107 -0.02261879 0.33566905
#> 512 -0.27375107 -0.06747884 0.42577436
#> 513 -0.25333326 -0.06747884 0.46887864
#> 514 -0.25333326 -0.09701684 0.53656227
#> 515 -0.25333326 -0.03517485 0.39429922
#> 516 -0.25333326 -0.08625276 0.51417389
#> 517 -0.25333326 -0.04822813 0.42030071
#> 518 -0.25333326 -0.04966244 0.41524434
#> 519 -0.23572265 -0.04966244 0.45227102
#> 520 -0.23572265 -0.04966244 0.45227102
#> 521 -0.23117368 -0.04966244 0.45833590
#> 522 -0.22441138 -0.04966244 0.47060919
#> 523 -0.21960818 -0.04966244 0.47848849
#> 524 -0.14716006 -0.04966244 0.69424313
#> 525 -0.14716006 -0.04603695 0.67838086
#> 526 -0.16736889 -0.04603695 0.61634713
#> 527 -0.16736889 -0.06564889 0.67037369
#> 528 -0.16736889 -0.02859086 0.56039998
#> 529 -0.17267290 -0.02859086 0.54135254
#> 530 -0.17267290 -0.01660923 0.50237360
#> 531 -0.17267290 -0.04887628 0.59273030
#> 532 -0.17267290 -0.04949283 0.58934493
#> 533 -0.17267290 -0.02381784 0.51148186
#> 534 -0.17267290 -0.03138729 0.52806820
#> 535 -0.17267290 -0.02518923 0.50512824
#> 536 -0.17267290 -0.03497588 0.52933485
#> 537 -0.19412491 -0.03497588 0.46430319
#> 538 -0.19412491  0.01944371 0.33628952
#> 539 -0.19412491  0.03135457 0.30532484
#> 540 -0.20333675  0.03135457 0.28116702
#> 541 -0.20333675  0.04236117 0.25462777
#> 542 -0.20175774  0.04236117 0.25250298
#> 543 -0.20175774  0.06894365 0.20343285
#> 544 -0.18170348  0.06894365 0.23586429
#> 545 -0.21957411  0.06894365 0.17553605
#> 546 -0.24518562  0.06894365 0.13954054
#> 547 -0.24518562  0.09476043 0.10902193
#> 548 -0.24518562  0.09527378 0.10481956
#> 549 -0.24518562  0.08945477 0.10738478
#> 550 -0.25739017  0.08945477 0.09253380
#> 551 -0.25739017  0.08945477 0.09253380
#> 552 -0.25739017  0.08945477 0.09253380
#> 553 -0.25739017  0.08945477 0.09253380
#> 554 -0.25739017  0.08945477 0.09253380
#> 555 -0.25739017  0.08945477 0.09253380
#> 556 -0.25739017  0.08945477 0.09253380
#> 557 -0.25739017  0.08945477 0.09253380
#> 558 -0.25739017  0.08945477 0.09253380
#> 559 -0.25739017  0.08945477 0.09253380
#> 560 -0.25739017  0.08945477 0.09253380
#> 561 -0.25739017  0.08945477 0.09253380
#> 562 -0.25739017  0.08945477 0.09253380
#> 563 -0.25739017  0.08945477 0.09253380
#> 564 -0.25739017  0.08945477 0.09253380
#> 565 -0.25739017  0.08945477 0.09253380
#> 566 -0.25739017  0.08945477 0.09253380
#> 567 -0.25739017  0.08945477 0.09253380
#> 568 -0.25739017  0.08945477 0.09253380
#> 569 -0.25739017  0.08945477 0.09253380
#> 570 -0.25739017  0.08945477 0.09253380
#> 571 -0.25739017  0.08945477 0.09253380
#> 572 -0.25739017  0.08945477 0.09253380
#> 573 -0.25739017  0.08945477 0.09253380
#> 574 -0.25739017  0.08945477 0.09253380
#> 575 -0.25739017  0.08945477 0.09253380
#> 576 -0.25739017  0.08945477 0.09253380
#> 577 -0.25739017  0.08945477 0.09253380
#> 578 -0.25739017  0.08945477 0.09253380
#> 579 -0.25739017  0.08945477 0.09253380
#> 580 -0.25739017  0.08945477 0.09253380
#> 581 -0.25739017  0.08945477 0.09253380
#> 582 -0.25739017  0.08945477 0.09253380
#> 583 -0.25739017  0.08945477 0.09253380
#> 584 -0.25739017  0.08945477 0.09253380
#> 585 -0.25739017  0.08945477 0.09253380
#> 586 -0.25739017  0.08945477 0.09253380
#> 587 -0.25739017  0.08945477 0.09253380
#> 588 -0.25739017  0.08945477 0.09253380
#> 589 -0.25739017  0.08945477 0.09253380
#> 590 -0.25739017  0.08945477 0.09253380
#> 591 -0.25739017  0.08945477 0.09253380
#> 592 -0.25739017  0.08945477 0.09253380
#> 593 -0.25739017  0.08945477 0.09253380
#> 594 -0.25739017  0.08945477 0.09253380
#> 595 -0.25739017  0.08945477 0.09253380
#> 596 -0.25739017  0.08945477 0.09253380
#> 597 -0.25739017  0.08945477 0.09253380
#> 598 -0.25739017  0.08945477 0.09253380
#> 599 -0.25739017  0.08945477 0.09253380
#> 600 -0.25739017  0.08945477 0.09253380
#> 601 -0.25739017  0.08945477 0.09253380
#> 602 -0.25739017  0.08945477 0.09253380
#> 603 -0.25739017  0.08945477 0.09253380
#> 604 -0.25739017  0.08945477 0.09253380
#> 605 -0.25739017  0.08945477 0.09253380
#> 606 -0.25739017  0.08945477 0.09253380
#> 607 -0.25739017  0.08945477 0.09253380
#> 608 -0.25739017  0.08945477 0.09253380
#> 609 -0.25739017  0.08945477 0.09253380
#> 610 -0.25739017  0.08945477 0.09253380
#> 611 -0.25739017  0.08945477 0.09253380
#> 612 -0.25739017  0.08945477 0.09253380
#> 613 -0.25739017  0.08945477 0.09253380
#> 614 -0.25739017  0.08945477 0.09253380
#> 615 -0.25739017  0.08945477 0.09253380
#> 616 -0.25739017  0.08945477 0.09253380
#> 617 -0.25739017  0.08945477 0.09253380
#> 618 -0.25739017  0.08945477 0.09253380
#> 619 -0.25739017  0.08945477 0.09253380
#> 620 -0.25739017  0.08945477 0.09253380
#> 621 -0.25739017  0.08945477 0.09253380
#> 622 -0.25739017  0.08945477 0.09253380
#> 623 -0.25739017  0.08945477 0.09253380
#> 624 -0.25739017  0.08945477 0.09253380
#> 625 -0.25739017  0.08945477 0.09253380
#> 626 -0.25739017  0.08945477 0.09253380
#> 627 -0.25739017  0.08945477 0.09253380
#> 628 -0.25739017  0.08945477 0.09253380
#> 629 -0.25739017  0.08945477 0.09253380
#> 630 -0.25739017  0.08945477 0.09253380
#> 631 -0.25739017  0.08945477 0.09253380
#> 632 -0.25739017  0.08945477 0.09253380
#> 633 -0.25739017  0.08945477 0.09253380
#> 634 -0.25739017  0.08945477 0.09253380
#> 635 -0.25739017  0.08945477 0.09253380
#> 636 -0.25739017  0.08945477 0.09253380
#> 637 -0.25739017  0.08945477 0.09253380
#> 638 -0.25739017  0.08945477 0.09253380
#> 639 -0.25739017  0.08945477 0.09253380
#> 640 -0.25739017  0.08945477 0.09253380
#> 641 -0.25739017  0.08945477 0.09253380
#> 642 -0.25739017  0.08945477 0.09253380
#> 643 -0.25739017  0.08945477 0.09253380
#> 644 -0.25739017  0.08945477 0.09253380
#> 645 -0.25739017  0.08945477 0.09253380
#> 646 -0.25739017  0.08945477 0.09253380
#> 647 -0.25739017  0.08945477 0.09253380
#> 648 -0.25739017  0.08945477 0.09253380
#> 649 -0.25739017  0.08945477 0.09253380
#> 650 -0.25739017  0.08945477 0.09253380
#> 651 -0.25739017  0.08945477 0.09253380
#> 652 -0.25739017  0.08945477 0.09253380
#> 653 -0.25739017  0.08945477 0.09253380
#> 654 -0.25739017  0.08945477 0.09253380
#> 655 -0.25739017  0.08945477 0.09253380
#> 656 -0.25739017  0.08945477 0.09253380
#> 657 -0.25739017  0.08945477 0.09253380
#> 658 -0.25739017  0.08945477 0.09253380
#> 659 -0.25739017  0.08945477 0.09253380
#> 660 -0.25739017  0.08945477 0.09253380
#> 661 -0.25739017  0.08945477 0.09253380
#> 662 -0.25739017  0.08945477 0.09253380
#> 663 -0.25739017  0.08945477 0.09253380
#> 664 -0.25739017  0.08945477 0.09253380
#> 665 -0.25739017  0.08945477 0.09253380
#> 666 -0.25739017  0.08945477 0.09253380
#> 667 -0.25739017  0.08945477 0.09253380
#> 668 -0.25739017  0.08945477 0.09253380
#> 669 -0.25739017  0.08945477 0.09253380
#> 670 -0.25739017  0.08945477 0.09253380
#> 671 -0.25739017  0.08945477 0.09253380
#> 672 -0.25739017  0.08945477 0.09253380
#> 673 -0.25739017  0.08945477 0.09253380
#> 674 -0.25739017  0.08945477 0.09253380
#> 675 -0.25739017  0.08945477 0.09253380
#> 676 -0.25739017  0.08945477 0.09253380
#> 677 -0.25739017  0.08945477 0.09253380
#> 678 -0.25739017  0.08945477 0.09253380
#> 679 -0.25739017  0.08945477 0.09253380
#> 680 -0.25739017  0.08945477 0.09253380
#> 681 -0.25739017  0.08945477 0.09253380
#> 682 -0.25739017  0.08945477 0.09253380
#> 683 -0.25739017  0.08945477 0.09253380
#> 684 -0.25739017  0.08945477 0.09253380
#> 685 -0.25739017  0.08945477 0.09253380
#> 686 -0.25739017  0.08945477 0.09253380
#> 687 -0.25739017  0.08945477 0.09253380
#> 688 -0.25739017  0.08945477 0.09253380
#> 689 -0.25739017  0.08945477 0.09253380
#> 690 -0.25739017  0.08945477 0.09253380
#> 691 -0.25739017  0.08945477 0.09253380
#> 692 -0.25739017  0.08945477 0.09253380
#> 693 -0.25739017  0.08945477 0.09253380
#> 694 -0.25739017  0.08945477 0.09253380
#> 695 -0.25739017  0.08945477 0.09253380
#> 696 -0.25739017  0.08945477 0.09253380
#> 697 -0.25739017  0.08945477 0.09253380
#> 698 -0.25739017  0.08945477 0.09253380
#> 699 -0.25739017  0.08945477 0.09253380
#> 700 -0.25739017  0.08945477 0.09253380
#> 701 -0.25739017  0.08945477 0.09253380
#> 702 -0.25739017  0.08945477 0.09253380
#> 703 -0.25739017  0.08945477 0.09253380
#> 704 -0.25739017  0.08945477 0.09253380
#> 705 -0.25739017  0.08945477 0.09253380
#> 706 -0.25739017  0.08945477 0.09253380
#> 707 -0.25739017  0.08945477 0.09253380
#> 708 -0.25739017  0.08945477 0.09253380
#> 709 -0.25739017  0.08945477 0.09253380
#> 710 -0.25739017  0.08945477 0.09253380
#> 711 -0.25739017  0.08945477 0.09253380
#> 712 -0.25739017  0.08945477 0.09253380
#> 713 -0.25739017  0.08945477 0.09253380
#> 714 -0.25739017  0.08945477 0.09253380
#> 715 -0.25739017  0.08945477 0.09253380
#> 716 -0.25739017  0.08945477 0.09253380
#> 717 -0.25739017  0.08945477 0.09253380
#> 718 -0.25739017  0.08945477 0.09253380
#> 719 -0.25739017  0.08945477 0.09253380
#> 720 -0.25739017  0.08945477 0.09253380
#> 721 -0.25739017  0.08945477 0.09253380
#> 722 -0.25739017  0.08945477 0.09253380
#> 723 -0.25739017  0.08945477 0.09253380
#> 724 -0.25739017  0.08945477 0.09253380
#> 725 -0.25739017  0.08945477 0.09253380
#> 726 -0.25739017  0.08945477 0.09253380
#> 727 -0.25739017  0.08945477 0.09253380
#> 728 -0.25739017  0.08945477 0.09253380
#> 729 -0.25739017  0.08945477 0.09253380
#> 730 -0.25739017  0.08945477 0.09253380
#> 731 -0.25739017  0.08945477 0.09253380
#> 732 -0.25739017  0.08945477 0.09253380
#> 733 -0.25739017  0.08945477 0.09253380
#> 734 -0.25739017  0.08945477 0.09253380
#> 735 -0.25739017  0.08945477 0.09253380
#> 736 -0.25739017  0.08945477 0.09253380
#> 737 -0.25739017  0.08945477 0.09253380
#> 738 -0.25739017  0.08945477 0.09253380
#> 739 -0.25739017  0.08945477 0.09253380
#> 740 -0.25739017  0.08945477 0.09253380
#> 741 -0.25739017  0.08945477 0.09253380
#> 742 -0.25739017  0.08945477 0.09253380
#> 743 -0.25739017  0.08945477 0.09253380
#> 744 -0.25739017  0.08945477 0.09253380
#> 745 -0.25739017  0.08945477 0.09253380
#> 746 -0.25739017  0.08945477 0.09253380
#> 747 -0.25739017  0.08945477 0.09253380
```

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
prepends `replicate`, `timepoint`, and `analysis` columns to your result
data. Each row shows either `"interim"` or `"final"` in the `analysis`
column, reflecting which named analysis produced it. This naming
convention becomes essential in later examples where multiple named
analyses fire per replicate.
