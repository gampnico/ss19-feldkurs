
# Preliminary derivation of H under free convection conditions i.e
# little wind shear, strong instability
k = 0.4  # von Karman constant
cp = 1004  # J kg^-1 K^-1, heat capacity of air
r_d = 287.05  # J kg^-1 K^-1, specific gas constant for dry air
beta_1 = 0.86
g = 9.81

station = "schiessstand"
filename = "2019-05-24"

z_eff = pw.return_z_effective(station)
computed_data = free_flux_compute(filename)

print("Iteration started...\n")
for index, row in computed_data.iterrows():
    # Initialise first iteration step
    ob_ini = -1000
    l_ob = obukhov_iteration(index, ob_ini)
    computed_data.loc[index, "l_ob"] = l_ob
    step = 0
    while abs(l_ob - ob_ini) > 0.1:
        step += 1
        ob_prev = ob_ini
        ob_ini = l_ob
        l_ob = obukhov_iteration(index, ob_ini)
        # avoid infinite loops
        if math.isclose(l_ob, ob_prev, rel_tol=0.1) or step > 1000:
            break

print("Iteration completed!")
# Calculate sensible heat flux
computed_data["H"] = -computed_data["rho_air"] * cp \
                     * computed_data["theta_star"] * computed_data["u_star"]
# Calculate momentum flux
computed_data["momentum_flux"] = computed_data["rho_air"] * computed_data[
    "u_star"] ** 2
