library(here)
library(tidyverse)

harm1 <- read_csv(here("data/mae_reconstruction.csv"))
harm2 <- read_csv(here("data/mae_reconstruction_two_harmonics.csv"))
harm3 <- read_csv(here("data/mae_reconstruction_three_harmonics.csv"))
smooth <- read_csv(here("data/mae_reconstruction_smoothed.csv"))
interpolated <- read_csv(here("data/mae_reconstruction_interpolated.csv"))

dt_long <- harm1 |>
  pivot_longer(cols = everything()) |>
  mutate(harmonics = "One Harmonic") |>
  bind_rows(
    harm2 |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Two Harmonics")
  ) |>
  bind_rows(
    harm3 |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Three Harmonics")
  ) |>
  bind_rows(
    smooth |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Smoothed")
  ) |>
  bind_rows(
    interpolated |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Interpolated")
  )

dt_long |>
  filter(!str_detect(name, "original")) |>
  mutate(
    across(
      harmonics,
      ~factor(
        .x,
        levels = c(
          "One Harmonic", "Two Harmonics", "Three Harmonics",
          "Smoothed", "Interpolated"
          ),
        ordered = TRUE
        )
      )
    ) |>
  ggplot(aes(x = harmonics, y = value)) +
  geom_boxplot(fill = "dodgerblue1", alpha = 0.2) +
  labs(x = "", y = "MAE") +
  theme_bw()

quantile_df <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1)) {
  tibble(
    value = quantile(x, probs, na.rm = TRUE),
    prob = probs
  )
}

dt_long |>
  group_by(harmonics) |>
  reframe(across(value, quantile_df)) |>
  unnest_wider(value) |>
  pivot_wider(names_from = harmonics, values_from = value)

dt_long |>
  group_by(harmonics) |>
  summarise(
    across(value, list("mean" = mean, "sd" = sd))
    )

sharm1 <- read_csv(here("data/signal_reconstruction.csv"))
sharm2 <- read_csv(here("data/signal_reconstruction_two_harmonics.csv"))
sharm3 <- read_csv(here("data/signal_reconstruction_three_harmonics.csv"))
ssmooth <- read_csv(here("data/signal_reconstruction_smoothed.csv"))
sinterpolated <- read_csv(here("data/signal_reconstruction_interpolated.csv"))
truth <- read_csv(here("data/co2_signals.csv"))

index <- seq(from = 1, to = 144001, by = 10)

truth |>
  slice(index) |>
  mutate(t = seq(from = 0, by = 1/30, length.out = 14401)) |>
  pivot_longer(cols = -t) |>
  filter(name == "0009_8min.mat") |>
  mutate(signal = "truth") |>
  mutate(across(value, ~.x - mean(.x))) |>
  bind_rows(
    sinterpolated |>
      mutate(t = seq(from = 0, by = 1/30, length.out = 14401)) |>
      pivot_longer(cols = -t) |>
      filter(name == "0009_8min_INTERPOLATED.mat") |>
      mutate(signal = "recon")
  ) |>
  ggplot(aes(x = t, y = value)) +
  geom_line(aes(color = signal))

truth |>
  mutate(across(everything(), ~.x - mean(.x))) |>
  slice(index) |>
  select("0009_8min.mat") |>
  bind_cols(
    sharm1 |>
      select("0009_8min.mat")
  ) |>
  magrittr::set_names(c("truth", "recon")) |>
  ggplot(aes(x = recon, y = truth)) +
  geom_point(alpha = 0.3) +
  labs(
    y = "CO2",
    x = "Reconstruction",
    title = "Reconstruction Scatterplot",
    subtitle = "One Harmonic"
    ) +
  theme_bw()

truth |>
  mutate(across(everything(), ~.x - mean(.x))) |>
  slice(index) |>
  select("0009_8min.mat") |>
  bind_cols(
    sharm1 |>
      select("0009_8min.mat")
  ) |>
  magrittr::set_names(c("truth", "recon")) |>
  cor()

dt_long |>
  filter(harmonics == "One Harmonic") |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 13, fill = "dodgerblue4") +
  labs(x = "MAE", y = "", title = "MAE Distribution", subtitle = "One Harmonic") +
  theme_bw()

wide_harm1 <- read_csv(here("data/mae_reconstruction_wide.csv"))
wide_sharm1 <- read_csv(here("data/signal_reconstruction_wide.csv"))
wide_harm2 <- read_csv(here("data/mae_reconstruction_two_harmonics_wide.csv"))
wide_sharm2 <- read_csv(here("data/signal_reconstruction_two_harmonics_wide.csv"))

truth |>
  mutate(across(everything(), ~.x - mean(.x))) |>
  slice(index) |>
  select("0009_8min.mat") |>
  bind_cols(
    wide_sharm1 |>
      select("0009_8min.mat")
  ) |>
  magrittr::set_names(c("truth", "recon")) |>
  ggplot(aes(x = recon, y = truth)) +
  geom_point(alpha = 0.3) +
  labs(
    y = "CO2",
    x = "Reconstruction",
    title = "Reconstruction Scatterplot",
    subtitle = "One Harmonic"
  ) +
  theme_bw()

dt_long <- harm1 |>
  pivot_longer(cols = everything()) |>
  mutate(harmonics = "One Harmonic") |>
  bind_rows(
    harm2 |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Two Harmonics")
  ) |>
  bind_rows(
    harm3 |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Three Harmonics")
  ) |>
  bind_rows(
    smooth |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Smoothed")
  ) |>
  bind_rows(
    interpolated |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Interpolated")
  ) |>
  bind_rows(
    wide_harm1 |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "One/Wide")
  ) |>
  bind_rows(
    wide_harm2 |>
      pivot_longer(cols = everything()) |>
      mutate(harmonics = "Two/Wide")
  )

dt_long |>
  filter(!str_detect(name, "original")) |>
  mutate(
    across(
      harmonics,
      ~factor(
        .x,
        levels = c(
          "One Harmonic", "Two Harmonics", "Three Harmonics",
          "Smoothed", "Interpolated", "One/Wide", "Two/Wide"
        ),
        ordered = TRUE
      )
    )
  ) |>
  ggplot(aes(x = harmonics, y = value)) +
  geom_boxplot(fill = "dodgerblue1", alpha = 0.2) +
  labs(x = "", y = "MAE") +
  theme_bw()
