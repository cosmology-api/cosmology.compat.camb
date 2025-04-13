# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


def parse_package_authors(author_email):
    """Get names from a package's Author-email field."""
    import email, email.policy

    msg = email.message_from_string(f"To: {author_email}", policy=email.policy.default)
    return ", ".join(address.display_name for address in msg["to"].addresses)


# -- Project information -----------------------------------------------------

import importlib.metadata

metadata = importlib.metadata.metadata("cosmology.compat.camb")

project = metadata["Name"]
author = parse_package_authors(metadata["Author-email"])
copyright = f"2025, {author}"
release = metadata["Version"]
version = release.partition("-")[0]


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
# html_static_path = ["_static"]
