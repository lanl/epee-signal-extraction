#####VARIABLE SETTING#######
match_date1 <- as.POSIXct("2023-09-28 22:01:22", tz = "UTC")
match_date2 <- as.POSIXct("2023-09-29 05:00:00", tz = "UTC")

bin_select   <- paste0("bin", seq(1, 100, 1))
density_scale <- 1e-12
current_scale <- 1e-9
charge        <- 1.60218e-19
area          <- .000169964
velocity      <- 7600
oxygen_mass   <- 2.6566962e-26

####clean and merge all data
raw <- rawepee_process(calendar_dates, rawdata_pathname, julian_days)
energy_full <- raw$energy
current_full <- raw$current

energy  <- raw$energy  %>% dplyr::filter(dplyr::between(timestamp, match_date1, match_date2))
current <- raw$current %>% dplyr::filter(dplyr::between(timestamp, match_date1, match_date2))

llae      <- process_llae(lla_pathname, match_date1, match_date2)
epee_max  <- process_max(current, energy)
epee_max_full <- process_max(current_full, energy_full)
sv        <- process_sv(calendar_dates, sv_pathname, match_date1, match_date2, oxygen_mass, velocity, charge)
fpmu      <- process_fpmu(fpmu_pathname, match_date1, match_date2)

full_df <- merge_data(llae, epee_max, sv, fpmu)

#####FIGURES BEGIN HERE#######

p <- current %>%
  tidyr::pivot_longer(
    cols = starts_with("V"),
    names_to = "energy_bin",
    values_to = "current"
  ) %>%
  dplyr::mutate(
    energy_bin = as.numeric(sub("V", "", energy_bin))
  ) %>%
ggplot2::ggplot(., aes(x = timestamp, y = energy_bin, fill = current)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  scale_x_datetime(date_breaks = "1 hour", date_labels = "%m-%d:%H") +
  labs(
    x = NULL,
    y = "energy bin",
    fill = "current (nA)"
  ) +
  themes +
  ggplot2::theme(
    legend.position   = "bottom",
    axis.text.x       = ggplot2::element_text(angle = 0, hjust = 0.5),
    legend.key.width  = grid::unit(2, "cm"),
    legend.text       = ggplot2::element_text(size = 10)
  )
p
ggplot2::ggsave(here::here("figures", "paper", "fig_image.png"), plot = p, width = 12, height = 9, dpi = 300)


####FIGURE 1
x <- 1000
timedate_x <- current %>% dplyr::filter(timeind == x) %>% dplyr::pull(timestamp)
max_current <- apply(current[x, paste0("V", 1:100)], 1, max)
column_name <- colnames(current)[which.max(current[x, paste0("V", 1:100)])]
max_energy  <- energy[x, ] %>% dplyr::select(column_name) %>% dplyr::pull()

p <- energy[x, ] %>%
  tidyr::pivot_longer(., cols = c(V1:V100), names_to = "energy_bin", values_to = "e_value") %>%
  dplyr::full_join(
    (current[x, ] %>%
       tidyr::pivot_longer(., cols = c(V1:V100), names_to = "energy_bin", values_to = "c_value")),
    by = dplyr::join_by("timestamp", "energy_bin")
  ) %>%
  ggplot2::ggplot(data = ., ggplot2::aes(x = e_value, y = c_value)) +
  ggplot2::geom_point(ggplot2::aes(x = e_value, y = c_value, color = e_value), size = 3) +
  ggplot2::geom_line(color = "black", linewidth = 1) +
  ggplot2::geom_hline(yintercept = max_current, color = "grey", linetype = "dashed", linewidth = 1) +
  ggplot2::geom_vline(xintercept = max_energy,  color = "grey", linetype = "dashed", linewidth = 1) +
  ggplot2::scale_x_continuous(breaks = seq(0, 200, 50), limits = c(0, 200)) +
  ggplot2::scale_y_continuous(breaks = seq(0.0, 0.5, 0.1), limits = c(0, 0.5)) +
  ggplot2::labs(x = "energy (eV)", y = "current (nA)") +
  themes +
  ggplot2::theme(legend.position = "None") +
  ggplot2::geom_segment(
    ggplot2::aes(x = max_energy, y = max_current, xend = max_energy + 20, yend = max_current - 0.1),
    color = "darkblue", linetype = "dashed", linewidth = 1
  ) +
  ggplot2::geom_label(
    ggplot2::aes(x = max_energy + 20, y = max_current - 0.1, label = "italic(E['0,t'])~and~italic(I['0,t'])"),
    parse = TRUE, hjust = 0, size = 6, color = "darkblue", fill = "white", label.size = 1
  )

p
ggplot2::ggsave(here::here("figures", "paper", "fig1.png"), plot = p, width = 12, height = 9, dpi = 300)

##FIGURE 2
time_value_x <- 41556
timedate_x   <- current %>% dplyr::filter(timeind == time_value_x) %>% dplyr::pull(timestamp)
max_current  <- epee_max_full %>% dplyr::filter(utc_time_epee %in% timedate_x) %>% dplyr::pull(current_value_org)
column_name  <- epee_max_full %>% dplyr::filter(utc_time_epee %in% timedate_x) %>% dplyr::pull(energy_bin_org)
max_energy   <- epee_max_full %>% dplyr::filter(utc_time_epee %in% timedate_x) %>% dplyr::pull(energy_value_org)

p <- energy %>% dplyr::filter(timeind == time_value_x) %>%
  tidyr::pivot_longer(., cols = c(V1:V100), names_to = "energy_bin", values_to = "e_value") %>%
  dplyr::full_join(
    (current %>% dplyr::filter(timeind == time_value_x) %>%
       tidyr::pivot_longer(., cols = c(V1:V100), names_to = "energy_bin", values_to = "c_value")),
    by = dplyr::join_by("timestamp", "energy_bin")
  ) %>%
  ggplot2::ggplot(data = ., ggplot2::aes(x = e_value, y = c_value)) +
  ggplot2::geom_point(ggplot2::aes(x = e_value, y = c_value, color = e_value), size = 3) +
  ggplot2::geom_line(color = "black", linewidth = 1) +
  ggplot2::geom_hline(yintercept = max_current, color = "grey", linetype = "dashed", linewidth = 1) +
  ggplot2::geom_vline(xintercept = max_energy,  color = "grey", linetype = "dashed", linewidth = 1) +
  ggplot2::geom_segment(
    ggplot2::aes(x = max_energy - 50, y = max_current + 0.02, xend = max_energy, yend = max_current),
    color = "darkblue", linetype = "dashed", linewidth = 1
  ) +
  ggplot2::geom_label(
    ggplot2::aes(x = max_energy - 50, y = max_current + 0.02, label = "italic(E['0,t'])~and~italic(I['0,t'])"),
    parse = TRUE, hjust = 0, size = 6, color = "darkblue", fill = "white", label.size = 1
  ) +
  ggplot2::scale_x_continuous(breaks = seq(0, 200, 50), limits = c(0, 200)) +
  ggplot2::scale_y_continuous(breaks = seq(0.0, 0.2, 0.05), limits = c(0, 0.2)) +
  ggplot2::labs(x = "energy (eV)", y = "current (nA)") +
  themes +
  ggplot2::theme(legend.position = "None")

p
ggplot2::ggsave(here::here("figures", "paper", "fig2.png"), plot = p, width = 12, height = 9, dpi = 300)

##FIGURE 3
current_target <- 0.15

p <- epee_max %>%
  ggplot2::ggplot(., ggplot2::aes(x = utc_time_epee)) +
  ggplot2::geom_point(
    data = epee_max, ggplot2::aes(y = energy_value_org, color = current_value_org),
    size = 1, alpha = 0.5
  ) +
  ggplot2::geom_point(
    data = epee_max %>% dplyr::filter(current_value_org <= current_target),
    ggplot2::aes(y = energy_value_org),
    color = "darkgrey", size = 2
  ) +
  ggplot2::scale_color_gradient(low = "blue", high = "yellow", name = "current (nA)") +
  ggplot2::scale_x_datetime(date_labels = "%m-%d:%H", date_breaks = "1 hour") +
  ggplot2::scale_y_continuous(limits = c(0, 200), breaks = c(seq(0, 200, by = 50))) +
  ggplot2::labs(x = "", y = "energy (eV)") +
  themes +
  ggplot2::theme(
    legend.position   = "bottom",
    axis.text.x       = ggplot2::element_text(size = 10, angle = 0, vjust = 0.5),
    legend.key.width  = grid::unit(2, "cm"),
    legend.text       = ggplot2::element_text(size = 10)
  )
p
ggplot2::ggsave(here::here("figures", "paper", "fig3.png"), plot = p, width = 12, height = 9, dpi = 300)

##FIGURE 4 - process flow chart

#FIGURE 5
calendar_date_choice <- 1 #pertains to either Julian day 271 or 272
quarter_day <- 4
full_result <- estimate_noise_floor(current_full, energy_full, calendar_dates,
                                    calendar_date_choice = calendar_date_choice, quarter_day = quarter_day)


p1 <- ggplot2::ggplot(
  full_result$df_long %>% dplyr::filter(timeind %in% full_result$noisefloor_timeinds & energy_bin <= 97),
  ggplot2::aes(x = energy_bin, y = value, group = timeind, color = energy)
) +
  ggplot2::geom_point(size = 1, alpha = 0.5) +
  ggplot2::geom_line(
    data = data.frame(bin = 1:100, value = full_result$floor_profile) %>% dplyr::filter(bin <= 97),
    ggplot2::aes(x = bin, y = value),
    color = "red", inherit.aes = FALSE
  ) +
  ggplot2::scale_y_continuous(limits = c(0.02, 0.16), breaks = seq(0.02, 0.16, 0.02)) +
  ggplot2::labs(x = "", y = "current (nA)", color = "energy (eV)") +
  themes +
  ggplot2::scale_color_gradient(low = "#000000", high = "#007BFF", limits = c(0, 196), breaks = c(50, 100, 150), name = "energy (eV)") +
  ggplot2::theme(
    axis.text.x      = ggplot2::element_text(size = 10, angle = 0, vjust = 0.5),
    legend.position  = "bottom",
    legend.key.width = grid::unit(2, "cm"),
    legend.text      = ggplot2::element_text(size = 10),
    legend.title     = ggplot2::element_text(size = 12),
    legend.direction = "horizontal",
    legend.box       = "horizontal"
  )

quarter_day <- 1
calendar_date_choice <- 1
full_result <- estimate_noise_floor(current_full, energy_full, calendar_dates,
                                    calendar_date_choice = calendar_date_choice, quarter_day = quarter_day)

tmp <- full_result$df_long %>%
  dplyr::filter(timeind %in% full_result$noisefloor_timeinds) %>%
  dplyr::group_by(energy_bin) %>%
  dplyr::mutate(avg_predmu = mean(predmu, na.rm = TRUE)) %>%
  dplyr::ungroup()

p2 <- ggplot2::ggplot(
  full_result$df_long %>% dplyr::filter(timeind %in% full_result$noisefloor_timeinds & energy_bin <= 97),
  ggplot2::aes(x = energy_bin, y = value, group = timeind, color = energy)
) +
  ggplot2::geom_point(size = 1, alpha = 0.5) +
  ggplot2::geom_line(
    data = data.frame(bin = 1:100, value = full_result$floor_profile) %>% dplyr::filter(bin <= 97),
    ggplot2::aes(x = bin, y = value),
    color = "red", inherit.aes = FALSE
  ) +
  ggplot2::scale_y_continuous(limits = c(0.02, 0.16), breaks = seq(0.02, 0.16, 0.02)) +
  ggplot2::labs(x = "energy bin", y = " current (nA)", color = "energy (eV)") +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::scale_color_gradient(low = "#000000", high = "#007BFF", limits = c(0, 196), breaks = c(50, 100, 150), name = "energy (eV)") +
  themes +
  ggplot2::theme(
    legend.position  = "none",
    axis.text.x      = ggplot2::element_text(size = 10, angle = 0, vjust = 0.5),
    legend.key.width = grid::unit(2, "cm"),
    legend.text      = ggplot2::element_text(size = 10),
    legend.title     = ggplot2::element_text(size = 12),
    legend.direction = "horizontal",
    legend.box       = "horizontal"
  )

##combined p1/p2 plot
combined <- p1 / p2 + patchwork::plot_layout(guides = "collect") & ggplot2::theme(
  legend.position  = "bottom",
  legend.key.width = grid::unit(2, "cm"),
  legend.text      = ggplot2::element_text(size = 10),
  legend.title     = ggplot2::element_text(size = 12),
  legend.direction = "horizontal",
  legend.box       = "horizontal"
)
combined
ggplot2::ggsave(here::here("figures", "paper", "fig5.png"), plot = combined, width = 12, height = 9, dpi = 300)

######FIGURE 6
p <- full_df %>%
  ggplot2::ggplot(ggplot2::aes(x = utc_time_epee)) +
  ggplot2::geom_point(
    data = full_df %>% dplyr::filter(current_value_org < current_target),
    ggplot2::aes(x = utc_time_epee, y = energy_value_org),
    color = "darkgrey", size = 2
  ) +
  ggplot2::geom_point(
    data = full_df %>% dplyr::filter(!(is.na(energy_value_sv))),
    ggplot2::aes(y = energy_value_sv, color = current_value_sv),
    size = 1, alpha = 0.5
  ) +
  ggplot2::geom_point(
    data = full_df %>% dplyr::filter(!(is.na(energy_omitted))),
    ggplot2::aes(y = energy_value_org),
    color = "black", size = 3
  ) +
  ggplot2::scale_color_gradient(low = "blue", high = "yellow", name = "current (nA)") +
  ggplot2::scale_x_datetime(date_labels = "%m-%d:%H", date_breaks = "1 hour") +
  ggplot2::labs(x = "", y = "energy(eV)") +
  themes +
  ggplot2::theme(
    legend.position  = "bottom",
    axis.text.x      = ggplot2::element_text(size = 12, angle = 0, vjust = 0.5),
    legend.key.width = grid::unit(2, "cm"),
    legend.text      = ggplot2::element_text(size = 10)
  )
p
ggplot2::ggsave(here::here("figures", "paper", "fig6.png"), plot = p, width = 12, height = 9, dpi = 300)

#####FIGURE 7
p1 <- full_df %>%
  tidyr::pivot_longer(., cols = c(scp_epee_org, vf_fpp, vf_wlp), names_to = "process", values_to = "scp_value") %>%
  ggplot2::ggplot(ggplot2::aes(x = utc_time_epee, y = scp_value, color = process)) +
  ggplot2::geom_point(size = 1) +
  ggplot2::scale_x_datetime(date_labels = "%d:%H", date_breaks = "1 hour") +
  viridis::scale_color_viridis(
    discrete = TRUE, direction = -1, option = "D",
    labels = c(scp_epee_org = "original", vf_fpp = "FPMU FPP", vf_wlp = "FPMU WLP")
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 30), breaks = c(seq(0, 30, by = 5))) +
  ggplot2::labs(x = "", y = "spacecraft potential (V)", color = NULL) +
  themes +
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4))) +
  ggplot2::theme(
    legend.position  = "bottom",
    axis.text.x      = ggplot2::element_text(size = 12, angle = 0, vjust = 0.5),
    legend.key.width = grid::unit(2, "cm"),
    legend.text      = ggplot2::element_text(size = 10)
  )

p2 <- full_df %>%
  tidyr::pivot_longer(., cols = c(scp_epee_sv, vf_fpp, vf_wlp), names_to = "process", values_to = "scp_value") %>%
  ggplot2::ggplot(ggplot2::aes(x = utc_time_sv, y = scp_value, color = process)) +
  ggplot2::geom_point(size = 1) +
  ggplot2::scale_x_datetime(date_labels = "%d:%H", date_breaks = "1 hour") +
  viridis::scale_color_viridis(
    discrete = TRUE, direction = -1, option = "D",
    labels = c(scp_epee_sv = "smoothed", vf_fpp = "FPMU FPP", vf_wlp = "FPMU WLP")
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 30), breaks = c(seq(0, 30, by = 5))) +
  ggplot2::labs(x = "day:hour", y = "", color = NULL) +
  themes +
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4))) +
  ggplot2::theme(
    legend.position  = "bottom",
    axis.text.x      = ggplot2::element_text(size = 12, angle = 0, vjust = 0.5),
    legend.key.width = grid::unit(2, "cm"),
    legend.text      = ggplot2::element_text(size = 10)
  )

p3 <- p1 / p2
p3
ggplot2::ggsave(here::here("figures", "paper", "fig7.png"), plot = p3, width = 12, height = 9, dpi = 300)

#### FIGURE 8
p <- full_df %>%
  ggplot2::ggplot(ggplot2::aes(x = long, y = lat)) +
  ggplot2::geom_point(ggplot2::aes(color = scp_epee_sv), size = 2, alpha = 0.5) +
  ggplot2::geom_point(data = full_df %>% dplyr::filter(current_value_org < 0.15), color = "grey", size = 3) +
  ggplot2::geom_point(data = full_df %>% dplyr::filter(!is.na(current_omitted)), color = "black", size = 2) +
  ggplot2::scale_color_gradient(
    low = "blue", high = "yellow", name = "Spacecraft\nPotential (V)",
    limits = c(0, 30), breaks = seq(0, 30, by = 5)
  ) +
  ggplot2::scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 30)) +
  ggplot2::scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, 20)) +
  ggplot2::labs(x = "longitude", y = "latitude") +
  themes +
  ggplot2::theme(
    legend.position  = "bottom",
    axis.text.x      = ggplot2::element_text(size = 12, angle = 0, vjust = 0.5),
    legend.key.width = grid::unit(2, "cm"),
    legend.text      = ggplot2::element_text(size = 10),
    legend.title     = ggplot2::element_text(size = 12)
  )
p
ggplot2::ggsave(here::here("figures", "paper", "fig8.png"), plot = p, width = 12, height = 9, dpi = 300)

######END PAPER FIGURES######

###SUPPLEMENT FIGURES####

#######LOOKING AT QUARTER DAY 1 FOR DAY 2###########
quarter_day <- 1
calendar_date_choice <- 2
full_result <- estimate_noise_floor(
  current_full, energy_full, calendar_dates,
  calendar_date_choice = calendar_date_choice, quarter_day = quarter_day
)

tmp    <- full_result$tmp
tmp_bg <- full_result$tmp_bg

floor_profile     <- tmp$floor_profile
gaussian_profile  <- tmp$gaussian_profile
richards_profile  <- tmp$richards_profile
parabola_profile  <- tmp$parabola_profile
richards_params   <- tmp$richards_params
gaussian_params   <- tmp$gaussian_params
parabola_params   <- tmp$parabola_params

p <- data.frame(
  bin = 1:100,
  raw = tmp_bg,
  noise_floor = floor_profile,
  richards_profile = richards_profile,
  parabola_profile = parabola_profile,
  gaussian = gaussian_profile
) %>%
  ggplot2::ggplot(., ggplot2::aes(x = bin)) +
  ggplot2::geom_line(ggplot2::aes(y = c(tmp_bg[1:97], rep(NA, 3))), color = "steelblue1", linewidth = 1.3, linetype = "solid") +
  ggplot2::geom_line(ggplot2::aes(y = noise_floor), color = "black", linewidth = 2, linetype="dotted") +
  ggplot2::geom_line(ggplot2::aes(y = gaussian), color = "goldenrod", linetype = "dotdash", linewidth = 1.2) +
  ggplot2::geom_line(ggplot2::aes(y = richards_profile), color = "steelblue4", linewidth = 1.2) +
  ggplot2::geom_line(ggplot2::aes(y = parabola_profile), color = "magenta4", linewidth = 1.2) +
  ggplot2::scale_y_continuous(breaks = seq(0, 0.15, 0.05), limits = c(-0.01, 0.15)) +
  ggplot2::labs(x = "energy bin", y = "current (nA)") +
  themes +
  ggplot2::theme(
    axis.text.x   = ggplot2::element_text(size = 10, angle = 0, vjust = 0.5),
    plot.caption  = ggplot2::element_text(size = 10, hjust = 1, face = "italic"),
    legend.position = "none"
  )
p
ggplot2::ggsave(here::here("figures", "paper", "fig9.png"), plot = p, width = 12, height = 9, dpi = 300)


df <- data.frame(
  parameter = c(
      "Richards $A$", "Richards $k$", "Richards $x_{0}$", "Richards $\\nu$", "Richards $A_{0}$",
      "Gaussian amplitude", "Gaussian $\\mu$", "Gaussian $\\sigma$",
      "Parabola $p$", "Parabola $q$", "Parabola $r$"
    ),
  estimate  = round(c(richards_params, gaussian_params, parabola_params), 6)
) %>% dplyr::mutate(
  estimate = sprintf("%.6f", estimate),
  row = row_number(),
  highlight = case_when(
    row == 1        ~ "grey32",
    row == 6        ~ "goldenrod",
    row %in% 9:11   ~ "magenta4",
    row %in% 2:5    ~ "white",
    row %in% 7:8   ~ "white",
    TRUE            ~ NA_character_
  )
)

# long form for 2 columns: parameter (col 1, left), estimate (col 2, right)
tab_long <- df %>%
  dplyr::select(row, highlight, parameter, estimate) %>%
  tidyr::pivot_longer(cols = c(parameter, estimate), names_to = "col", values_to = "text") %>%
  dplyr::mutate(
    col_index = ifelse(col == "parameter", 1L, 2L),
    y = nrow(df) - row + 1L,
    # Convert TeX only for the left column; keep numbers as-is
    label = ifelse(
      col_index == 1L,
      vapply(
        text,
        function(s) paste(deparse(latex2exp::TeX(s)[[1]]), collapse = ""),
        character(1)
      ),
      text
    ),
    # text color rule: colored fill -> white text; white or NA fill -> black text
    font_col = dplyr::case_when(
      is.na(highlight)     ~ "black",
      highlight == "white" ~ "black",
      TRUE                 ~ "white"
    )
  )


n_rows <- nrow(df)

p_tbl <- ggplot2::ggplot(tab_long, aes(x = col_index, y = y)) +
  ggplot2::geom_tile(aes(fill = highlight),
                     width = 0.98, height = 0.98, color = "grey85", size = 0.3, show.legend = FALSE) +
  ggplot2::scale_fill_identity() +
  ggplot2::geom_text(
    aes(label = label,
        hjust = ifelse(col_index == 1L, 0, 1),
        fontface = ifelse(col_index == 1L, "bold", "plain"),
        colour = font_col),
    size = 3.2,
    nudge_x = ifelse(tab_long$col_index == 1L, -0.35, 0.35),
    parse = TRUE
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::annotate("rect", xmin = 0.5, xmax = 2.5, ymin = n_rows + 0.25, ymax = n_rows + 1.25,
                    fill = "grey94", color = "grey85", size = 0.3) +
  ggplot2::annotate("text", x = 1, y = n_rows + 0.75, label = "parameter",
                    fontface = "bold", hjust = 0, size = 3.2, nudge_x = -0.35) +
  ggplot2::annotate("text", x = 2, y = n_rows + 0.75, label = "estimate",
                    fontface = "bold", hjust = 1, size = 3.2, nudge_x =  0.35) +
  ggplot2::coord_cartesian(xlim = c(0.5, 2.5), ylim = c(0.5, n_rows + 1.25), clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(8, 12, 8, 12))

p_tbl
ggsave(here::here("figures", "paper", "fig10.png"), plot = p_tbl, width = 7, height = 6, dpi = 300)
