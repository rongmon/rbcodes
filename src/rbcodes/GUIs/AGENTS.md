# Version & Maintenance Notes
- Current version: v2.0.0
- Removed in v2.0.0: rb_plot_spec.py (superseded by multispecviewer)

# Commands
- Run all tests: `pytest`
- Run single test: `pytest tests/path/to/test.py::TestClass::test_function`
- Run with coverage: `pytest` (coverage report is default)
- Lint: `flake8` (Black-compatible config)
- Build docs: `tox -e docs`
- Run doctests: `tox -e doctests`

# Code Style
- Python >=3.9.6 required
- Black-compatible max line length: 88 chars
- Testing: pytest with unittest style classes
- Docstrings: Google style for all public functions/classes
- Imports order: stdlib, third-party (numpy etc), local
- Local imports use absolute paths (from rbcodes.GUIs.abstools import X)
- Naming: snake_case for functions/vars, PascalCase for classes
- Error handling: Explicit try/except with proper error types
- Core deps: PyQt5>=5.15.7, numpy>=1.22.3, linetools>=0.3
- Type hints: Optional but recommended for new code
- Use setuptools for packaging (setup.cfg for config)

# Project Structure 
- GUIs: zgui/, specgui/, multispecviewer/, abstools/
- Core logic in IGM/, halo/, utils/
- Tests in tests/ directory matching implementation paths
- Entry points: rb-multispec, rb-spec, rb-llsfitter, rb-zgui
- Documentation in docs/ with rst format