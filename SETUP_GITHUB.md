# Publishing this repository on GitHub

The repository is already initialised with two tagged commits (`v1.0.0`, `v1.1.0`).
Three things to do before the first push.

## 1. Replace the placeholders

`ORG` stands in for your GitHub owner (user or organisation) in the README badges,
`pyproject.toml` URLs, `CITATION.cff` and `CONTRIBUTING.md`:

```bash
grep -rl 'ORG/somascan-qc' . --exclude-dir=.git | xargs sed -i '' 's|ORG/somascan-qc|YOUR-ORG/somascan-qc|g'   # macOS
# grep -rl 'ORG/somascan-qc' . --exclude-dir=.git | xargs sed -i 's|ORG/somascan-qc|YOUR-ORG/somascan-qc|g'    # Linux
```

Then edit the author block in `CITATION.cff` (`GIVEN-NAME`, `FAMILY-NAME`,
`AFFILIATION`, optionally `orcid`) - GitHub renders a "Cite this repository" button
from that file, and it is the only file with deliberate placeholders left in it.

## 2. Set the commit identity

The two existing commits were made with a neutral identity. To claim them as your own:

```bash
git config user.name  "Your Name"
git config user.email "you@example.org"
git rebase --root --exec 'git commit --amend --no-edit --reset-author'
git tag -f -a v1.0.0 <first-commit-sha> -m "v1.0.0"
git tag -f -a v1.1.0 HEAD -m "v1.1.0"
```

If you would rather start from a clean single commit, delete `.git` and
`git init` again - nothing else in the repository depends on its history.

## 3. Create the remote and push

```bash
# with the GitHub CLI
gh repo create YOUR-ORG/somascan-qc --public --source=. --remote=origin \
   --description "Reproducible QC and filtering for SomaScan ADAT data"
git push -u origin main --follow-tags

# or manually, after creating an empty repo in the web UI
git remote add origin git@github.com:YOUR-ORG/somascan-qc.git
git push -u origin main --follow-tags
```

## 4. After the first push

- **Actions** run automatically on push and pull request: tests on Python 3.9-3.12
  (Linux) and 3.12 (macOS), lint, an end-to-end demo run, a byte-level reproducibility
  check and a Docker build that executes the demo inside the image. Check the badge
  turns green; nothing needs enabling.
- **Topics**: `proteomics`, `somascan`, `somalogic`, `quality-control`, `aptamer`,
  `reproducible-research`.
- **Release**: `gh release create v1.1.0 --notes-file CHANGELOG.md` (or push the tag and
  draft it in the UI).
- **Archiving for citation (optional)**: link the repository to Zenodo before cutting a
  release and Zenodo will mint a DOI for it; add the DOI badge to the README and the
  `doi:` field to `CITATION.cff`.
- **Branch protection (optional)**: require the `tests` and `docker build` checks on
  `main`.

## What is deliberately not in the repository

No ADAT files, no result matrices, no reports - `.gitignore` excludes `*.adat`,
`*.tsv`, `qc_report.html` and the output directories so study data cannot be committed
by accident. The test suite and CI run on a simulated ADAT generated at runtime by
`somascan_qc/simulate.py`, so the repository stays free of real study data.
