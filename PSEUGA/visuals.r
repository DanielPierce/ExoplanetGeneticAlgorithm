library("car")

run_name = "kepler8testing"

fitness_csv_path <- paste(run_name, "/", run_name, "_FITNESS.csv", sep = "", collapse = NULL)
best_curve_csv_path <- paste(run_name, "/", "best_curve.csv", sep = "", collapse = NULL)
target_curve_csv_path <- paste(run_name, "/", "target_curve.csv", sep = "", collapse = NULL)

fitnesses <- read.csv(fitness_csv_path, header = FALSE, col.names = c("XCorrMax (^)", "XCorrDistFromCenter (v)", "MinDistance (v)", "group"))
best_curve <- read.csv(best_curve_csv_path, header = FALSE, row.names = c("Seconds from epoch", "Flux (e/s)"))
target_curve <- read.csv(target_curve_csv_path, header = FALSE, row.names = c("Seconds from epoch", "Flux (e/s)"))

xcorr <- fitnesses$XCorrMax....
xcorr_dist <- fitnesses$XCorrDistFromCenter..v.
min_dist <- fitnesses$MinDistance..v.


scatter3d(x = xcorr, y = xcorr_dist, z = min_dist, surface = FALSE, groups = as.factor(fitnesses$group), grid = TRUE)

best_timings <- as.numeric(best_curve[1,])
best_fluxes <- as.numeric(best_curve[2,])

target_timings <- as.numeric(target_curve[1,])
target_fluxes <- as.numeric(target_curve[2,])

#plot(target_timings, target_fluxes, type = 'l', xlab = "Seconds from Epoch", ylab = "Generated Flux (e/s)", xlim = c(min(target_timings),max(target_timings)), ylim = c(min(target_fluxes), max(target_fluxes)), col = "blue")

#lines(best_timings, best_fluxes, type = 'l', xlab = "Seconds from Epoch", ylab = "Generated Flux (e/s)", col = "red")