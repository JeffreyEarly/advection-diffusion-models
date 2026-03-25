Documentation
==============

+ The `WebsiteDocumentation` folder contains the hand-authored markdown that is copied into the repository `docs/` folder for GitHub Pages publishing.
+ The `../tools/build_website_documentation.m` script builds the published website by copying the website source, generating the version-history page from `../CHANGELOG.md`, and extracting class documentation with `class-docs`.
+ The generated `../docs` folder is committed output. Edit `Documentation/WebsiteDocumentation` and rerun the build script rather than editing `docs/` by hand.
