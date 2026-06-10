.PHONY: build run stop logs test lint clean

VENV  = .venv/bin
PORT ?= 8501

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
	$(VENV)/pytest tests/ -q

lint:
	$(VENV)/ruff check .

clean:
	docker compose down --rmi local -v --remove-orphans
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
