.PHONY: help install dev test lint demo docker docker-run clean
VERSION := 1.1.0
IMAGE   := somascan-qc:$(VERSION)

help:           ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

install:        ## Install the package into the active environment
	pip install .

dev:            ## Editable install with test and lint extras
	pip install -e ".[test,dev]"

test:           ## Run the test suite
	pytest -q

lint:           ## Static checks
	ruff check somascan_qc tests

demo:           ## Simulate an ADAT and run the pipeline on it
	somascan-qc demo -o somascan_qc_demo

docker:         ## Build the container image
	docker build -t $(IMAGE) .

docker-run:     ## Verify the image by running the demo inside it
	docker run --rm -v "$$PWD/somascan_qc_demo":/data $(IMAGE) demo -o /data

clean:
	rm -rf build dist *.egg-info .pytest_cache .ruff_cache somascan_qc_demo
	find . -name __pycache__ -type d -exec rm -rf {} +
