dump:
	@cd scripts && python dump_samples_with_aas.py IN

tsm:
	@cd scripts && python calc_chi_squares.py IN && python filter_holm_method.py IN
	@cd scripts && python calc_chi_squares.py IN B && python filter_holm_method.py IN B
