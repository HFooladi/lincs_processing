dose = [x.split() for x in query_trt.pert_idose]
time = [x.split() for x in query_trt.pert_itime]

pert_dose = list(map(lambda x: float(x[0]), dose))
pert_dose_unit = list(map(lambda x: x[1], dose))
pert_time = list(map(lambda x: float(x[0]), time))
pert_time_unit = list(map(lambda x: x[1], time))

adding = pd.DataFrame({
    'pert_dose': pert_dose,
    'pert_dose_unit': pert_dose_unit,
    'pert_time': pert_time,
    'pert_time_unit': pert_time_unit
})

adding = adding.set_index(query_trt.index)
query_trt = pd.concat([query_trt, adding], axis=1)
