"""Your list of terms: the terms that you suspect, which the report sets against the residuals, each against its twin,
in the question matched. It is empty at the start, and the package names none: with it empty the report reads the
residuals alone. A term is written here as a term file writes it, a term of the package that `python -m twinfit terms`
lists or a Mechanism of your own, and the names of model.py are read here as a term file reads them."""
from twinfit import mechanisms  # noqa: F401

TERMS = []
