#load in package for cleaning data
library(tidyverse)
library(fixest)
library(splines)
#package to create graphs
library(plotly)

#read in nthe 4 CSV files of NNDSS data, 8 pathogens in 4 files
nndss_test      <- read_csv("data/NNDSS_Cyclo_Campy.csv")
nndss_salmon    <- read_csv("data/NNDSS_Salmon_Listy.csv")
nndss_stec_shig <- read_csv("data/NNDSS_STEC_Shig.csv")
nndss_vib_yers  <- read_csv("data/NNDSS_Vib_Yers.csv")

#combine the datasets
nndss_all <- bind_rows(nndss_test, nndss_salmon, nndss_stec_shig, nndss_vib_yers)

valid_states <- c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island",
                  "Vermont", "New Jersey", "New York", "New York City", "Pennsylvania",
                  "Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin", "Iowa", "Kansas",
                  "Minnesota", "Missouri", "Nebraska", "North Dakota", "South Dakota",
                  "Delaware", "District of Columbia", "Florida", "Georgia", "Maryland",
                  "North Carolina", "South Carolina", "Virginia", "West Virginia",
                  "Alabama", "Kentucky", "Mississippi", "Tennessee", "Arkansas", "Louisiana",
                  "Oklahoma", "Texas", "Arizona", "Colorado", "Idaho", "Montana", "Nevada",
                  "New Mexico", "Utah", "Wyoming", "Alaska", "California", "Hawaii",
                  "Oregon", "Washington")


nndss_model_data <- nndss_all %>%
  
#turn the reporting area into "state" to standardize
  mutate(state = str_to_title(`Reporting Area`)) %>%
  filter(state %in% valid_states) %>%
  
#combine NYC and New York
  mutate(state = if_else(state == "New York City", "New York", state)) %>%
  mutate(Label = case_when(

#combine probable and confirmed listeriosis cases
    Label %in% c("Listeriosis, Confirmed", "Listeriosis, Probable") ~ "Listeriosis",
    str_detect(Label, "^Vibriosis") ~ "Vibriosis",
    str_detect(Label, "^Salmonellosis") ~ "Salmonellosis",
    Label == "Shiga toxin-producing Escherichia coli (STEC)" ~ "STEC",
    TRUE ~ Label
  )) %>%
  mutate(week_date = as.Date(paste0(`Current MMWR Year`, "-W", sprintf("%02d", `MMWR WEEK`), "-1"),
                             format = "%Y-W%U-%u")) %>%
  filter(!is.na(week_date)) %>%
  group_by(state, Label, week_date) %>%
  summarise(cases = sum(`Current week`, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    downgraded = if_else(Label %in% c("Campylobacteriosis", "Cyclosporiasis", "Listeriosis",
                                      "Shigellosis", "Vibriosis"), 1, 0),
    post_policy = if_else(week_date >= as.Date("2025-07-01"), 1, 0),
    week_of_year = as.numeric(format(week_date, "%U"))
  )

nndss_model_data %>% count(Label)   

nndss_national <- nndss_model_data %>%
  group_by(Label, week_date) %>%
  summarise(total_cases = sum(cases), .groups = "drop")

# standardize dates from CSV files to Year, Month, Date format
ggplot(nndss_national, aes(x = week_date, y = total_cases, color = Label)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2025-07-01"), linetype = "dashed", color = "black") +
  labs(title = "Weekly Reported Cases by Pathogen",
       subtitle = "Dashed line = July 1, 2025 policy change",
       x = "Week", y = "Total reported cases (all states)") +
  theme_minimal()

did_model_full <- feols(cases ~ post_policy * downgraded | state + week_date,
                        data = nndss_model_data,
                        cluster = "state")
summary(did_model_full)

run_pathogen_did <- function(pathogen_name, data) {
  model_data <- data %>%
    filter(Label %in% c(pathogen_name, "Salmonellosis", "STEC")) %>%
    mutate(downgraded = if_else(Label == pathogen_name, 1, 0))
  
  feols(cases ~ post_policy * downgraded | state + week_date,
        data = model_data,
        cluster = "state")
}

pathogens <- c("Campylobacteriosis", "Cyclosporiasis", "Listeriosis", "Shigellosis", "Vibriosis")

results <- map(pathogens, ~ run_pathogen_did(.x, nndss_model_data))
names(results) <- pathogens

walk(names(results), ~ {
  cat("\n\n====", .x, "====\n")
  print(summary(results[[.x]]))
})


estimate_undercount <- function(pathogen_name, data) {
  path_data <- data %>%
    filter(Label == pathogen_name) %>%
#To answer: how many cases were reported in x state during # week
    group_by(week_date, week_of_year) %>%
    summarise(cases = sum(cases), .groups = "drop")
  
  pre_data <- path_data %>% filter(week_date < as.Date("2025-07-01"))
  trend_model <- lm(cases ~ ns(week_of_year, df = 8) + week_date, data = pre_data)
  
  path_data %>%
    mutate(predicted_cases = predict(trend_model, newdata = path_data))
}

undercount_results <- map_dfr(pathogens, ~ estimate_undercount(.x, nndss_model_data) %>%
                                mutate(pathogen = .x))

#calculates the gap between expected and actual cases
undercount_summary <- undercount_results %>%
  filter(week_date >= as.Date("2025-07-01")) %>%
  group_by(pathogen) %>%
  summarise(
    actual_cases = sum(cases),
    expected_cases = sum(predicted_cases),
    gap = expected_cases - actual_cases,
    pct_undercounted = round(gap / expected_cases * 100, 1)
  )

print(undercount_summary, width = Inf)

ggplot(undercount_results, aes(x = week_date)) +
  geom_line(aes(y = cases, color = "Actual reported")) +
  geom_line(aes(y = predicted_cases, color = "Expected (pre-policy trend)"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2025-07-01"), linetype = "dotted", color = "black") +
  facet_wrap(~ pathogen, scales = "free_y") +
  labs(title = "Actual vs. Expected Cases by Pathogen",
       subtitle = "Dashed line = counterfactual trend if policy hadn't changed",
       x = "Week", y = "Cases", color = "") +
  theme_minimal()


estimate_undercount_by_state <- function(pathogen_name, data) {
  path_data <- data %>% filter(Label == pathogen_name)
  
  pre_data_national <- path_data %>%
    filter(week_date < as.Date("2025-07-01")) %>%
    group_by(week_date, week_of_year) %>%
    summarise(cases = sum(cases), .groups = "drop")
  
  trend_model <- lm(cases ~ ns(week_of_year, df = 8) + week_date, data = pre_data_national)
  
  state_share <- path_data %>%
    filter(week_date < as.Date("2025-07-01")) %>%
    group_by(state) %>%
    summarise(state_total = sum(cases)) %>%
    mutate(share = state_total / sum(state_total))
  
  post_dates <- path_data %>%
    filter(week_date >= as.Date("2025-07-01")) %>%
    distinct(week_date, week_of_year)
  
  total_expected_national <- sum(predict(trend_model, newdata = post_dates))
  
  path_data %>%
    filter(week_date >= as.Date("2025-07-01")) %>%
    group_by(state) %>%
    summarise(actual_cases = sum(cases), .groups = "drop") %>%
    left_join(state_share, by = "state") %>%
    mutate(
      expected_cases = share * total_expected_national,
      gap = expected_cases - actual_cases,
      pct_undercounted = round(gap / expected_cases * 100, 1),
      pathogen = pathogen_name
    )
}

state_undercount <- map_dfr(pathogens, ~ estimate_undercount_by_state(.x, nndss_model_data))

state_undercount_clean <- state_undercount %>%
  filter(is.finite(pct_undercounted), state_total >= 5) %>%
  mutate(state_abb = state.abb[match(state, state.name)])



plot_list <- lapply(pathogens, function(p) {
  df <- state_undercount_clean %>% filter(pathogen == p)
  list(
    type = 'choropleth',
    locationmode = 'USA-states',
    locations = df$state_abb,
    z = df$pct_undercounted,
    text = paste0(df$state, "<br>Actual: ", round(df$actual_cases),
                  "<br>Expected: ", round(df$expected_cases),
                  "<br>Undercount: ", df$pct_undercounted, "%"),
    hoverinfo = 'text',
    colorscale = 'RdBu',
    reversescale = TRUE,
    zmid = 0,
    visible = (p == pathogens[1])
  )
})

fig <- plot_ly()
for (tr in plot_list) fig <- add_trace(fig, type = tr$type, locationmode = tr$locationmode,
                                       locations = tr$locations, z = tr$z, text = tr$text,
                                       hoverinfo = tr$hoverinfo, colorscale = tr$colorscale,
                                       reversescale = tr$reversescale, zmid = tr$zmid,
                                       visible = tr$visible)

buttons <- lapply(seq_along(pathogens), function(i) {
  visible <- rep(FALSE, length(pathogens))
  visible[i] <- TRUE
  list(method = "restyle", args = list("visible", visible), label = pathogens[i])
})


fig <- fig %>% plotly::layout(
  title = list(
    text = "Did the CDC's Optional Pathogen Reporting Policy Lead to Fewer Detected Cases?",
    font = list(size = 18)
  ),
  geo = list(scope = 'usa'),
  updatemenus = list(list(y = 1, x = 0.1, buttons = buttons)),
  margin = list(t = 100, b = 180, l = 60, r = 60),
  annotations = list(
    list(
      text = "This map estimates how far reported cases of five pathogens fell based on pre-policy<br>trends since July 2025, when the CDC announced that it would make tracking of foodborne<br>illnesses from these pathogens optional. Some states are excluded due to unreliable low case counts.",
      xref = "paper", yref = "paper",
      x = 0.5, y = -0.15,
      xanchor = "center", yanchor = "top",
      showarrow = FALSE,
      font = list(size = 11, color = "gray40"),
      align = "center"
    )
  )
)

fig

library(htmlwidgets)
saveWidget(fig, "index.html", selfcontained = TRUE)