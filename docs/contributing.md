#### Updating mkdocs

First, install the development dependencies (includes mkdocs):

```bash
uv sync
```

Then use mkdocs via `uv run` (this ensures the project venv is used):

- `uv run mkdocs build` — build the site locally and check for errors
- `uv run mkdocs serve` — serve the site locally at http://127.0.0.1:8000 for live preview

Pushing to `main` (or merging a PR) automatically deploys the docs site via
the `Deploy Docs` GitHub Actions workflow, which also generates the pdoc API
reference. There's no need to deploy manually.
