function test_find_max_xy(evt)
	qmax, iqmax = findmax(evt.charge)
	simax = evt.sensor_id[iqmax]
	qm, sim = find_max_xy(evt,"sensor_id", "charge")
	return qm ≈ qmax && sim ≈ simax
end

function test_find_xyz_sipm_qmax(hitdf)
	qxmax, xmax = find_max_xy(hitdf,"x", "q")
	qymax, ymax = find_max_xy(hitdf,"y", "q")
	qzmax, zmax = find_max_xy(hitdf,"z", "q")
	@test qxmax ≈ qymax
	@test qxmax ≈ qzmax
end
