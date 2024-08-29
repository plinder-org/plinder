# Release

Code is contributed via pull requests on *GitHub*.
A new `plinder` version is released on each pull request merge.
In consequence each merge automatically induces a
[Semantic Version](https://semver.org/lang/de/) bump.
The version bumping semantics is controlled via the commit history since the previous
release:

- If `bumpversion skip` is present in the commit message, the version will not be bumped
- If `bumpversion major` is present in the commit message, the major version will be bumped
- If `bumpversion minor` is present in the commit message, the minor version will be bumped
- If `bumpversion patch` is present in the commit message (or nothing is found), the patch version will be bumped

:::{note}
The CI workflow will use the **most recent** match in the commit history to make its decision.
:::

Each new version release automatically triggers the following platforms:

- A new *PyPI* release is created.
- A new *Docker* image is pushed to the
  [registry](https://github.com/plinder-org/plinder/pkgs/container/plinder).
- This documentation website is updated.
