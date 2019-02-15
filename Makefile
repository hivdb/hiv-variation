dump:
	@cd scripts && python dump_samples_with_aas.py IN

in_tsm:
	@cd scripts && python calc_chi_squares.py IN && python filter_holm_method.py IN
	@cd scripts && python calc_chi_squares.py IN A && python filter_holm_method.py IN A
	@cd scripts && python calc_chi_squares.py IN B && python filter_holm_method.py IN B
	@cd scripts && python calc_chi_squares.py IN C && python filter_holm_method.py IN C
	@cd scripts && python calc_chi_squares.py IN CRF01_AE && python filter_holm_method.py IN CRF01_AE
	@cd scripts && python calc_chi_squares.py IN CRF02_AG && python filter_holm_method.py IN CRF02_AG

in_prevalence:
	DATABASE_URI="mysql+pymysql://rshafer:rshafer@10.77.6.244/HIVDB2" python scripts/export_aapcnt.py INSTI local/in_prevalence-no_filter.json
	DATABASE_URI="mysql+pymysql://rshafer:rshafer@10.77.6.244/HIVDB2" python scripts/export_aapcnt.py INSTI local/in_prevalence-no_qa_issues.json --filter NO_QA_ISSUES
	#@cd scripts && python calc_prevalence_by_subtypes.py IN
