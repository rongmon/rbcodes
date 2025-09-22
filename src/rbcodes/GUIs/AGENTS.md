# Commands
- Run all tests: `pytest`
- Run single test: `pytest tests/path/to/test.py::TestClass::test_function`  
- Run with coverage: `pytest` (coverage report is enabled by default)
- Lint: `flake8` (Black-compatible config, max line length 88)
- Build docs: `cd docs && make html`

# Code Style
- Python >=3.9.6 required
- Follow Black code style (max line length: 88 chars)
- Use absolute imports (from rbcodes.GUIs.abstools import X)
- Import order: stdlib, third-party (numpy etc), local 
- Docstrings: Google style for all public functions/classes
- Naming: snake_case for functions/vars, PascalCase for classes
- Error handling: Explicit try/except with proper error types
- Type hints recommended for new code
- Testing: pytest with unittest style classes
- Key deps: PyQt5>=5.15.7, numpy>=1.22.3, linetools>=0.3