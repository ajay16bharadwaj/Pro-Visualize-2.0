.PHONY: build run stop logs test lint clean

VENV  = .venv/bin
PORT ?= 8501

# Prefer the local venv binaries when present (developer machine); fall back to
# PATH so the same targets work in CI / Docker where there is no .venv.
PYTEST := $(if $(wildcard $(VENV)/pytest),$(VENV)/pytest,pytest)
RUFF   := $(if $(wildcard $(VENV)/ruff),$(VENV)/ruff,ruff)

build:
	docker compose build

run: build
	docker compose up -d
	@echo "App running at http://localhost:$(PORT)"

stop:
	docker compose down

logs:
	docker compose logs -f

test:
	$(PYTEST) tests/ -q

lint:
	$(RUFF) check .

clean:
	docker compose down --rmi local -v --remove-orphans
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
