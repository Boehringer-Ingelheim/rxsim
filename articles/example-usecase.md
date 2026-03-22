# Example usecase

## Setup

``` r
library(rxsim)
```

## Population

``` r
ann <- Population$new(name = "Ann", vector_to_dataframe(rnorm(50)))
enroll <- cbind(1:10, rep(5L, 10))
drop <- cbind(1:10, rep(1L, 10))

for (i in 1:nrow(enroll)) {
  ann$set_enrolled(time = enroll[i, 1], n = enroll[i, 2])
  ann$set_dropped(time = drop[i, 1], n = drop[i, 2])
}
```

## Timer

``` r
# --- Usage ---
# Example reader funcs
summary_reader <- function(dat, t_now) {
  list(
    t = t_now,
    n = nrow(dat),
    mean_value = if ("value" %in% names(dat) &&
      nrow(dat) > 0) {
      mean(dat$value)
    } else {
      NA_real_
    }
  )
}
ids_reader <- function(dat, t_now) {
  unique(dat$id)
}

# Sample data
set.seed(1)
df <- data.frame(
  id      = rep(1:5, each = 4),
  visit   = rep(1:4, times = 5),
  status  = sample(c("active", "inactive"), 20, replace = TRUE),
  value   = rnorm(20),
  center  = sample(c("EU", "US"), 20, replace = TRUE)
)

timers <- Timer$new("TidyTimers")
# Each reader carries its own dplyr-like predicates
timers$add_condition(status == "active",
  visit >= 3,
  analysis = summary_reader,
  name = "active_v3_sum"
)
timers$add_condition(center == "EU", value > 0, analysis = ids_reader, name = "eu_ids")

out <- timers$check_conditions(df, Sys.time())
```

## Trial

``` r
set.seed(123)

data <- rbind(
  cbind(
    id = 1:20,
    readout_time = 1,
    value = rnorm(20, mean = 50)
  ),
  cbind(
    id = 1:20,
    readout_time = 4,
    value = rnorm(20, mean = 50)
  )
)
pop1 <- Population$new("Arm A", data = data.frame(data))

timepoints <- data.frame(
  time = c(1,2,3.1,4,5,6),
  arm = rep("Arm A", 6),
  dropper = c(2L, rep(1L, 5)),
  enroller = rep(3L, 6)
)

# --- Timers with multiple timepoints ---
t <- Timer$new(name = "TrialTimers") # Use your updated Timers from earlier
add_timepoints(t, timepoints)

# t$add_condition(
#   length(value) > 4,
#   analysis = function(d, tt) {
#     mean(d$value)
#   },
#   name = "overall_mean"
# )

t$add_condition(
  sum((time - enroll_time) > 1) > 3,
  analysis = NULL,
  name = "study_time_>_2_for_3"
)

# t$add_condition(
#   (time - enroll_time) > 1,
#   analysis = NULL,
#   name = "blah"
# )

trial <- Trial$new(
  name = "Trial A",
  timer = t,
  population = list(pop1)
)

# --- Run ---
trial$run()
#> Warning:  returning filtered data as is because condition 'study_time_>_2_for_3' has no applicable analysis
#> Warning:  returning filtered data as is because condition 'study_time_>_2_for_3' has no applicable analysis
#> Warning:  returning filtered data as is because condition 'study_time_>_2_for_3' has no applicable analysis
#> Warning:  returning filtered data as is because condition 'study_time_>_2_for_3' has no applicable analysis

prettify_results(trial$results)
#>     time study_time_._2_for_3.id study_time_._2_for_3.readout_time
#> 1    3.1                       1                                 1
#> 2    3.1                       2                                 1
#> 3    3.1                       3                                 1
#> 4    3.1                       6                                 1
#> 5    3.1                       8                                 1
#> 6    3.1                       9                                 1
#> 7    3.1                      11                                 1
#> 8    3.1                      18                                 1
#> 9    3.1                      20                                 1
#> 10   3.1                       1                                 4
#> 11   3.1                       2                                 4
#> 12   3.1                       3                                 4
#> 13   3.1                       6                                 4
#> 14   3.1                       8                                 4
#> 15   3.1                       9                                 4
#> 16   3.1                      11                                 4
#> 17   3.1                      18                                 4
#> 18   3.1                      20                                 4
#> 19   4.0                       1                                 1
#> 20   4.0                       2                                 1
#> 21   4.0                       3                                 1
#> 22   4.0                       6                                 1
#> 23   4.0                       7                                 1
#> 24   4.0                       8                                 1
#> 25   4.0                       9                                 1
#> 26   4.0                      11                                 1
#> 27   4.0                      12                                 1
#> 28   4.0                      13                                 1
#> 29   4.0                      18                                 1
#> 30   4.0                      20                                 1
#> 31   4.0                       1                                 4
#> 32   4.0                       2                                 4
#> 33   4.0                       3                                 4
#> 34   4.0                       6                                 4
#> 35   4.0                       7                                 4
#> 36   4.0                       8                                 4
#> 37   4.0                       9                                 4
#> 38   4.0                      11                                 4
#> 39   4.0                      12                                 4
#> 40   4.0                      13                                 4
#> 41   4.0                      18                                 4
#> 42   4.0                      20                                 4
#> 43   5.0                       1                                 1
#> 44   5.0                       2                                 1
#> 45   5.0                       3                                 1
#> 46   5.0                       4                                 1
#> 47   5.0                       6                                 1
#> 48   5.0                       7                                 1
#> 49   5.0                       8                                 1
#> 50   5.0                       9                                 1
#> 51   5.0                      11                                 1
#> 52   5.0                      12                                 1
#> 53   5.0                      13                                 1
#> 54   5.0                      14                                 1
#> 55   5.0                      16                                 1
#> 56   5.0                      18                                 1
#> 57   5.0                      20                                 1
#> 58   5.0                       1                                 4
#> 59   5.0                       2                                 4
#> 60   5.0                       3                                 4
#> 61   5.0                       4                                 4
#> 62   5.0                       6                                 4
#> 63   5.0                       7                                 4
#> 64   5.0                       8                                 4
#> 65   5.0                       9                                 4
#> 66   5.0                      11                                 4
#> 67   5.0                      12                                 4
#> 68   5.0                      13                                 4
#> 69   5.0                      14                                 4
#> 70   5.0                      16                                 4
#> 71   5.0                      18                                 4
#> 72   5.0                      20                                 4
#> 73   6.0                       1                                 1
#> 74   6.0                       2                                 1
#> 75   6.0                       3                                 1
#> 76   6.0                       4                                 1
#> 77   6.0                       6                                 1
#> 78   6.0                       7                                 1
#> 79   6.0                       8                                 1
#> 80   6.0                       9                                 1
#> 81   6.0                      10                                 1
#> 82   6.0                      11                                 1
#> 83   6.0                      12                                 1
#> 84   6.0                      13                                 1
#> 85   6.0                      14                                 1
#> 86   6.0                      15                                 1
#> 87   6.0                      16                                 1
#> 88   6.0                      17                                 1
#> 89   6.0                      18                                 1
#> 90   6.0                      20                                 1
#> 91   6.0                       1                                 4
#> 92   6.0                       2                                 4
#> 93   6.0                       3                                 4
#> 94   6.0                       4                                 4
#> 95   6.0                       6                                 4
#> 96   6.0                       7                                 4
#> 97   6.0                       8                                 4
#> 98   6.0                       9                                 4
#> 99   6.0                      10                                 4
#> 100  6.0                      11                                 4
#> 101  6.0                      12                                 4
#> 102  6.0                      13                                 4
#> 103  6.0                      14                                 4
#> 104  6.0                      15                                 4
#> 105  6.0                      16                                 4
#> 106  6.0                      17                                 4
#> 107  6.0                      18                                 4
#> 108  6.0                      20                                 4
#>     study_time_._2_for_3.value study_time_._2_for_3.arm
#> 1                     49.43952                    Arm A
#> 2                     49.76982                    Arm A
#> 3                     51.55871                    Arm A
#> 4                     51.71506                    Arm A
#> 5                     48.73494                    Arm A
#> 6                     49.31315                    Arm A
#> 7                     51.22408                    Arm A
#> 8                     48.03338                    Arm A
#> 9                     49.52721                    Arm A
#> 10                    48.93218                    Arm A
#> 11                    49.78203                    Arm A
#> 12                    48.97400                    Arm A
#> 13                    48.31331                    Arm A
#> 14                    50.15337                    Arm A
#> 15                    48.86186                    Arm A
#> 16                    50.42646                    Arm A
#> 17                    49.93809                    Arm A
#> 18                    49.61953                    Arm A
#> 19                    49.43952                    Arm A
#> 20                    49.76982                    Arm A
#> 21                    51.55871                    Arm A
#> 22                    51.71506                    Arm A
#> 23                    50.46092                    Arm A
#> 24                    48.73494                    Arm A
#> 25                    49.31315                    Arm A
#> 26                    51.22408                    Arm A
#> 27                    50.35981                    Arm A
#> 28                    50.40077                    Arm A
#> 29                    48.03338                    Arm A
#> 30                    49.52721                    Arm A
#> 31                    48.93218                    Arm A
#> 32                    49.78203                    Arm A
#> 33                    48.97400                    Arm A
#> 34                    48.31331                    Arm A
#> 35                    50.83779                    Arm A
#> 36                    50.15337                    Arm A
#> 37                    48.86186                    Arm A
#> 38                    50.42646                    Arm A
#> 39                    49.70493                    Arm A
#> 40                    50.89513                    Arm A
#> 41                    49.93809                    Arm A
#> 42                    49.61953                    Arm A
#> 43                    49.43952                    Arm A
#> 44                    49.76982                    Arm A
#> 45                    51.55871                    Arm A
#> 46                    50.07051                    Arm A
#> 47                    51.71506                    Arm A
#> 48                    50.46092                    Arm A
#> 49                    48.73494                    Arm A
#> 50                    49.31315                    Arm A
#> 51                    51.22408                    Arm A
#> 52                    50.35981                    Arm A
#> 53                    50.40077                    Arm A
#> 54                    50.11068                    Arm A
#> 55                    51.78691                    Arm A
#> 56                    48.03338                    Arm A
#> 57                    49.52721                    Arm A
#> 58                    48.93218                    Arm A
#> 59                    49.78203                    Arm A
#> 60                    48.97400                    Arm A
#> 61                    49.27111                    Arm A
#> 62                    48.31331                    Arm A
#> 63                    50.83779                    Arm A
#> 64                    50.15337                    Arm A
#> 65                    48.86186                    Arm A
#> 66                    50.42646                    Arm A
#> 67                    49.70493                    Arm A
#> 68                    50.89513                    Arm A
#> 69                    50.87813                    Arm A
#> 70                    50.68864                    Arm A
#> 71                    49.93809                    Arm A
#> 72                    49.61953                    Arm A
#> 73                    49.43952                    Arm A
#> 74                    49.76982                    Arm A
#> 75                    51.55871                    Arm A
#> 76                    50.07051                    Arm A
#> 77                    51.71506                    Arm A
#> 78                    50.46092                    Arm A
#> 79                    48.73494                    Arm A
#> 80                    49.31315                    Arm A
#> 81                    49.55434                    Arm A
#> 82                    51.22408                    Arm A
#> 83                    50.35981                    Arm A
#> 84                    50.40077                    Arm A
#> 85                    50.11068                    Arm A
#> 86                    49.44416                    Arm A
#> 87                    51.78691                    Arm A
#> 88                    50.49785                    Arm A
#> 89                    48.03338                    Arm A
#> 90                    49.52721                    Arm A
#> 91                    48.93218                    Arm A
#> 92                    49.78203                    Arm A
#> 93                    48.97400                    Arm A
#> 94                    49.27111                    Arm A
#> 95                    48.31331                    Arm A
#> 96                    50.83779                    Arm A
#> 97                    50.15337                    Arm A
#> 98                    48.86186                    Arm A
#> 99                    51.25381                    Arm A
#> 100                   50.42646                    Arm A
#> 101                   49.70493                    Arm A
#> 102                   50.89513                    Arm A
#> 103                   50.87813                    Arm A
#> 104                   50.82158                    Arm A
#> 105                   50.68864                    Arm A
#> 106                   50.55392                    Arm A
#> 107                   49.93809                    Arm A
#> 108                   49.61953                    Arm A
#>     study_time_._2_for_3.enroll_time study_time_._2_for_3.drop_time
#> 1                                2.0                             NA
#> 2                                3.1                             NA
#> 3                                3.1                             NA
#> 4                                1.0                            2.0
#> 5                                1.0                            1.0
#> 6                                2.0                            3.1
#> 7                                1.0                            1.0
#> 8                                2.0                             NA
#> 9                                3.1                             NA
#> 10                               2.0                             NA
#> 11                               3.1                             NA
#> 12                               3.1                             NA
#> 13                               1.0                            2.0
#> 14                               1.0                            1.0
#> 15                               2.0                            3.1
#> 16                               1.0                            1.0
#> 17                               2.0                             NA
#> 18                               3.1                             NA
#> 19                               2.0                            4.0
#> 20                               3.1                             NA
#> 21                               3.1                             NA
#> 22                               1.0                            2.0
#> 23                               4.0                             NA
#> 24                               1.0                            1.0
#> 25                               2.0                            3.1
#> 26                               1.0                            1.0
#> 27                               4.0                             NA
#> 28                               4.0                             NA
#> 29                               2.0                             NA
#> 30                               3.1                             NA
#> 31                               2.0                            4.0
#> 32                               3.1                             NA
#> 33                               3.1                             NA
#> 34                               1.0                            2.0
#> 35                               4.0                             NA
#> 36                               1.0                            1.0
#> 37                               2.0                            3.1
#> 38                               1.0                            1.0
#> 39                               4.0                             NA
#> 40                               4.0                             NA
#> 41                               2.0                             NA
#> 42                               3.1                             NA
#> 43                               2.0                            4.0
#> 44                               3.1                             NA
#> 45                               3.1                             NA
#> 46                               5.0                             NA
#> 47                               1.0                            2.0
#> 48                               4.0                             NA
#> 49                               1.0                            1.0
#> 50                               2.0                            3.1
#> 51                               1.0                            1.0
#> 52                               4.0                             NA
#> 53                               4.0                             NA
#> 54                               5.0                             NA
#> 55                               5.0                             NA
#> 56                               2.0                            5.0
#> 57                               3.1                             NA
#> 58                               2.0                            4.0
#> 59                               3.1                             NA
#> 60                               3.1                             NA
#> 61                               5.0                             NA
#> 62                               1.0                            2.0
#> 63                               4.0                             NA
#> 64                               1.0                            1.0
#> 65                               2.0                            3.1
#> 66                               1.0                            1.0
#> 67                               4.0                             NA
#> 68                               4.0                             NA
#> 69                               5.0                             NA
#> 70                               5.0                             NA
#> 71                               2.0                            5.0
#> 72                               3.1                             NA
#> 73                               2.0                            4.0
#> 74                               3.1                             NA
#> 75                               3.1                             NA
#> 76                               5.0                             NA
#> 77                               1.0                            2.0
#> 78                               4.0                             NA
#> 79                               1.0                            1.0
#> 80                               2.0                            3.1
#> 81                               6.0                             NA
#> 82                               1.0                            1.0
#> 83                               4.0                             NA
#> 84                               4.0                             NA
#> 85                               5.0                             NA
#> 86                               6.0                             NA
#> 87                               5.0                             NA
#> 88                               6.0                             NA
#> 89                               2.0                            5.0
#> 90                               3.1                            6.0
#> 91                               2.0                            4.0
#> 92                               3.1                             NA
#> 93                               3.1                             NA
#> 94                               5.0                             NA
#> 95                               1.0                            2.0
#> 96                               4.0                             NA
#> 97                               1.0                            1.0
#> 98                               2.0                            3.1
#> 99                               6.0                             NA
#> 100                              1.0                            1.0
#> 101                              4.0                             NA
#> 102                              4.0                             NA
#> 103                              5.0                             NA
#> 104                              6.0                             NA
#> 105                              5.0                             NA
#> 106                              6.0                             NA
#> 107                              2.0                            5.0
#> 108                              3.1                            6.0
#>     study_time_._2_for_3.subject_id study_time_._2_for_3.measurement_time
#> 1                                 1                                   3.0
#> 2                                 2                                   4.1
#> 3                                 3                                   4.1
#> 4                                 4                                   2.0
#> 5                                 5                                   2.0
#> 6                                 6                                   3.0
#> 7                                 7                                   2.0
#> 8                                 8                                   3.0
#> 9                                 9                                   4.1
#> 10                                1                                   6.0
#> 11                                2                                   7.1
#> 12                                3                                   7.1
#> 13                                4                                   5.0
#> 14                                5                                   5.0
#> 15                                6                                   6.0
#> 16                                7                                   5.0
#> 17                                8                                   6.0
#> 18                                9                                   7.1
#> 19                                1                                   3.0
#> 20                                2                                   4.1
#> 21                                3                                   4.1
#> 22                                4                                   2.0
#> 23                                5                                   5.0
#> 24                                6                                   2.0
#> 25                                7                                   3.0
#> 26                                8                                   2.0
#> 27                                9                                   5.0
#> 28                               10                                   5.0
#> 29                               11                                   3.0
#> 30                               12                                   4.1
#> 31                                1                                   6.0
#> 32                                2                                   7.1
#> 33                                3                                   7.1
#> 34                                4                                   5.0
#> 35                                5                                   8.0
#> 36                                6                                   5.0
#> 37                                7                                   6.0
#> 38                                8                                   5.0
#> 39                                9                                   8.0
#> 40                               10                                   8.0
#> 41                               11                                   6.0
#> 42                               12                                   7.1
#> 43                                1                                   3.0
#> 44                                2                                   4.1
#> 45                                3                                   4.1
#> 46                                4                                   6.0
#> 47                                5                                   2.0
#> 48                                6                                   5.0
#> 49                                7                                   2.0
#> 50                                8                                   3.0
#> 51                                9                                   2.0
#> 52                               10                                   5.0
#> 53                               11                                   5.0
#> 54                               12                                   6.0
#> 55                               13                                   6.0
#> 56                               14                                   3.0
#> 57                               15                                   4.1
#> 58                                1                                   6.0
#> 59                                2                                   7.1
#> 60                                3                                   7.1
#> 61                                4                                   9.0
#> 62                                5                                   5.0
#> 63                                6                                   8.0
#> 64                                7                                   5.0
#> 65                                8                                   6.0
#> 66                                9                                   5.0
#> 67                               10                                   8.0
#> 68                               11                                   8.0
#> 69                               12                                   9.0
#> 70                               13                                   9.0
#> 71                               14                                   6.0
#> 72                               15                                   7.1
#> 73                                1                                   3.0
#> 74                                2                                   4.1
#> 75                                3                                   4.1
#> 76                                4                                   6.0
#> 77                                5                                   2.0
#> 78                                6                                   5.0
#> 79                                7                                   2.0
#> 80                                8                                   3.0
#> 81                                9                                   7.0
#> 82                               10                                   2.0
#> 83                               11                                   5.0
#> 84                               12                                   5.0
#> 85                               13                                   6.0
#> 86                               14                                   7.0
#> 87                               15                                   6.0
#> 88                               16                                   7.0
#> 89                               17                                   3.0
#> 90                               18                                   4.1
#> 91                                1                                   6.0
#> 92                                2                                   7.1
#> 93                                3                                   7.1
#> 94                                4                                   9.0
#> 95                                5                                   5.0
#> 96                                6                                   8.0
#> 97                                7                                   5.0
#> 98                                8                                   6.0
#> 99                                9                                  10.0
#> 100                              10                                   5.0
#> 101                              11                                   8.0
#> 102                              12                                   8.0
#> 103                              13                                   9.0
#> 104                              14                                  10.0
#> 105                              15                                   9.0
#> 106                              16                                  10.0
#> 107                              17                                   6.0
#> 108                              18                                   7.1
#>     study_time_._2_for_3.time
#> 1                         3.1
#> 2                         3.1
#> 3                         3.1
#> 4                         3.1
#> 5                         3.1
#> 6                         3.1
#> 7                         3.1
#> 8                         3.1
#> 9                         3.1
#> 10                        3.1
#> 11                        3.1
#> 12                        3.1
#> 13                        3.1
#> 14                        3.1
#> 15                        3.1
#> 16                        3.1
#> 17                        3.1
#> 18                        3.1
#> 19                        4.0
#> 20                        4.0
#> 21                        4.0
#> 22                        4.0
#> 23                        4.0
#> 24                        4.0
#> 25                        4.0
#> 26                        4.0
#> 27                        4.0
#> 28                        4.0
#> 29                        4.0
#> 30                        4.0
#> 31                        4.0
#> 32                        4.0
#> 33                        4.0
#> 34                        4.0
#> 35                        4.0
#> 36                        4.0
#> 37                        4.0
#> 38                        4.0
#> 39                        4.0
#> 40                        4.0
#> 41                        4.0
#> 42                        4.0
#> 43                        5.0
#> 44                        5.0
#> 45                        5.0
#> 46                        5.0
#> 47                        5.0
#> 48                        5.0
#> 49                        5.0
#> 50                        5.0
#> 51                        5.0
#> 52                        5.0
#> 53                        5.0
#> 54                        5.0
#> 55                        5.0
#> 56                        5.0
#> 57                        5.0
#> 58                        5.0
#> 59                        5.0
#> 60                        5.0
#> 61                        5.0
#> 62                        5.0
#> 63                        5.0
#> 64                        5.0
#> 65                        5.0
#> 66                        5.0
#> 67                        5.0
#> 68                        5.0
#> 69                        5.0
#> 70                        5.0
#> 71                        5.0
#> 72                        5.0
#> 73                        6.0
#> 74                        6.0
#> 75                        6.0
#> 76                        6.0
#> 77                        6.0
#> 78                        6.0
#> 79                        6.0
#> 80                        6.0
#> 81                        6.0
#> 82                        6.0
#> 83                        6.0
#> 84                        6.0
#> 85                        6.0
#> 86                        6.0
#> 87                        6.0
#> 88                        6.0
#> 89                        6.0
#> 90                        6.0
#> 91                        6.0
#> 92                        6.0
#> 93                        6.0
#> 94                        6.0
#> 95                        6.0
#> 96                        6.0
#> 97                        6.0
#> 98                        6.0
#> 99                        6.0
#> 100                       6.0
#> 101                       6.0
#> 102                       6.0
#> 103                       6.0
#> 104                       6.0
#> 105                       6.0
#> 106                       6.0
#> 107                       6.0
#> 108                       6.0
```

``` r

# --- Two populations with a common 'value' column ---
# long_format
pop1 <- Population$new("Arm A", data = data.frame(
  id = 1:20,
  value = rnorm(20, mean = 50),
  readout_time = 1
))
pop2 <- Population$new("Arm B", data = data.frame(
  id = 21:40,
  value = rnorm(20, mean = 55),
  readout_time = 1
))

# --- Timers with multiple timepoints ---
t <- Timer$new(name = "TrialTimers") # Use your updated Timers from earlier
add_timepoints(t, timepoints)

timepoints$time <- 1:6
timepoints$arm <- rep("Arm B", 6)
add_timepoints(t, timepoints)

# --- Reader conditions (per-reader predicates like dplyr::filter) ---
# IMPORTANT: func must accept (data, current_time). Wrap base functions accordingly.

# event condition
t$add_condition(
  length(value) > 4,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean"
)

t$add_condition(
  sum(value > 40) > 10,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean"
)

trigger_by_calendar(
  4, t,
  analysis = function(d, tt) {
    mean(d$value)
  }
)

t$add_condition(
  sum((value + time) > 45) > 10,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean_2"
)

# --- Trial with list of populations ---
trial <- Trial$new(
  name = "Trial A",
  timer = t,
  population = list(pop1, pop2)
)

# --- Run ---
trial$run()

prettify_results(trial$results)
#>   time overall_mean overall_mean_2 cal_time_4
#> 1  1.0     52.74784             NA         NA
#> 2  2.0     52.77916       52.77916         NA
#> 3  3.0     53.20100       53.20100         NA
#> 4  3.1     52.52427       52.52427         NA
#> 5  4.0     52.53041       52.53041   52.53041
#> 6  5.0     52.53819       52.53819         NA
#> 7  6.0     52.53907       52.53907         NA
```
