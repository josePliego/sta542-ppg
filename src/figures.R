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
  ggplot(aes(x = harmonics, y = value)) +
  geom_boxplot()

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
    across(value, list("mean" = mean, "sd" = sd)))
    )
