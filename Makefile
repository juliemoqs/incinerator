.PHONY: all pytest mypy ruff-check ruff-check-fix ruff-format isort

CMD:=uv run
PYMODULE:=src
TESTS:=tests

# Run all the checks which do not change files
all: pytest mypy ruff-check ruff-check-fix ruff-format build serve deploy

# Run tests using pytest
pytest:
	$(CMD) pytest $(PYMODULE) $(TESTS)

# Run tests using mypy
mypy:
	$(CMD) mypy $(PYMODULE)

# Ruff 
ruff-check:
	$(CMD) ruff check $(PYMODULE) $(TESTS)

# Ruff 
ruff-check-fix:
	$(CMD) ruff check --fix $(PYMODULE) $(TESTS)

# Ruff 
ruff-format:
	$(CMD) ruff format $(PYMODULE) $(TESTS)

# build docs
build:
	$(CMD) mkdocs build

serve:
	$(CMD) mkdocs serve

deploy:
	$(CMD) mkdocs gh-deploy --force