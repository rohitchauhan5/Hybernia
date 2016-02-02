for base_cbf = [45, 50, 55, 60, 65, 70, 75] % loop through CBF
    for base_cmr = [90, 110, 130, 150, 170, 190] % loop through CMR
        plot_all_give_pars(base_cbf, base_cmr) %input cbf in ml/100g/min, cmr in umol/100g/min
    end
end