local/in_prevalence.json: scripts/export_aapcnt.py
	@mkdir -p local/
	# Create local/in_prevalence.json
	DATABASE_URI="mysql+pymysql://rshafer:rshafer@10.77.6.244/HIVDB2" python scripts/export_aapcnt.py INSTI local/in_prevalence.json

local/in_prevalence-allow_mixture.json: scripts/export_aapcnt.py
	@mkdir -p local/
	# Create local/in_prevalence-allow_mixture.json
	DATABASE_URI="mysql+pymysql://rshafer:rshafer@10.77.6.244/HIVDB2" python scripts/export_aapcnt.py INSTI local/in_prevalence-allow_mixture.json --allow-mixture

local/in_drm_prevalence.txt: local/in_prevalence.json scripts/create_prevalence_table.py
	# Create local/in_drm_prevalence.txt
	@cd scripts && python create_prevalence_table.py IN -i ../local/in_prevalence.json -o ../local/in_drm_prevalence.txt

local/in_drm_prevalence-allow_mixture.txt: local/in_prevalence-allow_mixture.json scripts/create_prevalence_table.py
	# Create local/in_drm_prevalence-allow_mixture.txt
	@cd scripts && python create_prevalence_table.py IN -i ../local/in_prevalence-allow_mixture.json -o ../local/in_drm_prevalence-allow_mixture.txt

local/in_sigmut_prevalence.txt: local/in_prevalence.json scripts/create_prevalence_table.py
	# Create local/in_sigmut_prevalence.txt
	@cd scripts && python create_prevalence_table.py IN -i ../local/in_prevalence.json -o ../local/in_sigmut_prevalence.txt --num-algs-range 0 0

local/in_sigmut_prevalence-allow_mixture.txt: local/in_prevalence-allow_mixture.json scripts/create_prevalence_table.py
	# Create local/in_sigmut_prevalence-allow_mixture.txt
	@cd scripts && python create_prevalence_table.py IN -i ../local/in_prevalence-allow_mixture.json -o ../local/in_sigmut_prevalence-allow_mixture.txt --num-algs-range 0 0

local/in_chi2.txt: local/in_prevalence.json scripts/calc_chi_squares.py
	# Create local/in_chi_squares.txt
	@cd scripts && python calc_chi_squares.py IN -i ../local/in_prevalence.json -o ../local/in_chi2.txt

local/in_chi2-allow_mixture.txt: local/in_prevalence.json scripts/calc_chi_squares.py
	# Create local/in_chi_squares.txt
	@cd scripts && python calc_chi_squares.py IN -i ../local/in_prevalence-allow_mixture.json -o ../local/in_chi2-allow_mixture.txt

all: local/in_drm_prevalence.txt local/in_drm_prevalence-allow_mixture.txt local/in_sigmut_prevalence.txt local/in_sigmut_prevalence-allow_mixture.txt

.PHONY: all
