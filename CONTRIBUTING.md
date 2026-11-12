# Contributing to hacettepe-stem-cell

Thank you for your interest in contributing to **hacettepe-stem-cell**! We welcome contributions from biologists, data scientists, bioinformaticians, software engineers, and clinicians.

---

## 🛠️ Development Setup

1. **Fork and clone the repository:**
   ```bash
   git clone https://github.com/<your-username>/hacettepe-stem-cell.git
   cd hacettepe-stem-cell
   ```

2. **Set up a Python virtual environment:**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. **Install dependencies in editable mode:**
   ```bash
   pip install --upgrade pip
   pip install -r requirements.txt
   pip install -e .[dev]
   ```

4. **Run unit tests:**
   ```bash
   pytest -v
   ```

---

## 🧪 Testing Guidelines

* Every new algorithm, simulation engine, or data pipeline must include unit tests under `tests/`.
* We maintain >85% test coverage. Run coverage checks with:
  ```bash
  pytest --cov=src --cov-report=term-missing
  ```
* Ensure all tests pass across Python 3.10, 3.11, and 3.12.

---

## 📝 Code Style & Conventions

* **Format:** PEP 8 compliance. We use `black` for code formatting and `flake8` for linting.
* **Docstrings:** Use NumPy / Google style docstrings for all public classes, methods, and functions.
* **Typing:** Type hints are strongly encouraged (`typing.List`, `typing.Dict`, `typing.Optional`, etc.).
* **Commit Messages:** Follow conventional commits:
  * `feat: ...` for new features
  * `fix: ...` for bug fixes
  * `docs: ...` for documentation changes
  * `test: ...` for test additions or updates

---

## 🚀 Submitting a Pull Request

1. Create a feature branch:
   ```bash
   git checkout -b feature/my-new-algorithm
   ```
2. Commit your changes and push to your fork:
   ```bash
   git push origin feature/my-new-algorithm
   ```
3. Open a Pull Request against `main`.
4. Ensure CI checks pass. A maintainer will review your PR within a few days.
