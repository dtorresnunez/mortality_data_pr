#Esquema:
#Horario
#Modelo


formula_sb2 <- deaths ~
  factor(sex):period + #nuevo cambio: efecto de interacción sexo y período
  f(age_idx, model = model_age, constr = TRUE,
    hyper = list(theta = list(prior = SB2.prior(par_p_age, par_q_age, par_b_age)))) +
  f(region_idx, model = model_reg, graph = g, constr = TRUE,
    hyper = list(theta = list(prior = SB2.prior(par_p_reg , par_q_reg , par_b_reg)),
                 phi = list(prior = "logitbeta", param = c(0.5, 0.5)))) +
  f(period_idx, model = model_per, constr = TRUE,
    hyper = list(theta = list(prior = SB2.prior(par_p_per, par_q_per, par_b_per)))) +
  f(region_period_idx, model = model_s_t,
    hyper = list(theta = list(prior = SB2.prior(par_p_s_t, par_q_s_t, par_b_s_t)))) +
  f(cell_idx, model = model_cel,
    hyper = list(theta = list(prior = SB2.prior(par_p_cel, par_q_cel, par_b_cel))))
