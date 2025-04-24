# Makefile
env_name = mosaic-production

# Default target
.PHONY: help
help:
	@echo "Available commands:"
	@sed -n 's/^##//p' ${MAKEFILE_LIST} | column -t -s ':' | sed -e 's/^/ /'

## update_env: Update the conda environment
.PHONY: update_env
update_env:
	conda env update -n $(env_name) -f environment.mosaic-production-clis.yml --prune
	conda run -n $(env_name) python -m pip install -e ./earthdem-mosaic
	conda run -n $(env_name) python -m pip install -e ./rema-mosaic
	conda run -n $(env_name) python -m pip install -e ./matlib
	conda run -n $(env_name) python -m pip install -e ./mosaic-pipeline

## create_env: Create a new conda environment
.PHONY: create_env
create_env:
	conda env create -n $(env_name) -f environment.mosaic-production-clis.yml
	conda run -n $(env_name) python -m pip install -e ./earthdem-mosaic
	conda run -n $(env_name) python -m pip install -e ./rema-mosaic
	conda run -n $(env_name) python -m pip install -e ./matlib
	conda run -n $(env_name) python -m pip install -e ./mosaic-pipeline

## remove_env: Remove the conda environment
.PHONY: remove_env
remove_env:
	conda env remove -n $(env_name)

## activate_env: Show environment activation command
.PHONY: activate_env
activate_env:
	@echo "To activate the environment, use:"
	@echo ""
	@echo "conda activate $(env_name)"
	@echo ""