.PHONY: help install install-dev test lint format clean docker build

help:
	@echo "Phlox Development Commands"
	@echo "=========================="
	@echo "install        - Install package"
	@echo "install-dev    - Install package with dev dependencies"
	@echo "test           - Run tests with pytest"
	@echo "lint           - Run code linters"
	@echo "format         - Format code with black and isort"
	@echo "clean          - Remove build artifacts"
	@echo "docker         - Build Docker image"
	@echo "build          - Build distribution packages"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

test:
	pytest tests/ -v --cov=phlox --cov-report=html --cov-report=term

lint:
	flake8 src/phlox tests
	mypy src/phlox

format:
	black src/phlox tests
	isort src/phlox tests

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete

docker:
	docker build -t phlox:latest -f docker/Dockerfile .

build:
	python -m build
