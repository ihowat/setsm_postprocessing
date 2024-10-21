# Makefile
env_name = mosaic-production

# Default target
.PHONY: help
help:
	@echo "Available commands:"
	@sed -n 's/^##//p' ${MAKEFILE_LIST} | column -t -s ':' | sed -e 's/^/ /'

## update_env: Update the mamba environment
.PHONY: update_env
update_env:
	mamba env update -n $(env_name) -f environment.mosaic-production-clis.yml --prune
	mamba run -n $(env_name) python -m pip install -e ./earthdem-mosaic
	mamba run -n $(env_name) python -m pip install -e ./rema-mosaic

## create_env: Create a new mamba environment
.PHONY: create_env
create_env:
	mamba env create -n $(env_name) -f environment.mosaic-production-clis.yml
	mamba run -n $(env_name) python -m pip install -e ./earthdem-mosaic
	mamba run -n $(env_name) python -m pip install -e ./rema-mosaic

## remove_env: Remove the mamba environment
.PHONY: remove_env
remove_env:
	mamba env remove -n $(env_name)

## activate_env: Show environment activation command
.PHONY: activate_env
activate_env:
	@echo "To activate the environment, use:"
	@echo ""
	@echo "mamba activate $(env_name)"
	@echo ""