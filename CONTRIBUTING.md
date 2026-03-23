# Contributing to OpenQuantumSystem.jl

When contributing to this repository, please first discuss the change you wish
to make via issue, email, or any other method with the owners of this repository
before making a change. You can contribute with

- Reporting a bug
- Discussing the current state of the code and code review
- Submitting a fix
- Proposing new features

## Issues

We use GitHub issues to track public bugs. Please ensure your description is
clear and has sufficient instructions to be able to reproduce the issue.

## Setup

After cloning the repository, configure git to use the project's hooks:

```bash
git config core.hooksPath .githooks
```

This enables a `commit-msg` hook that blocks commits containing AI co-author attribution (e.g. `Co-Authored-By: Claude`).

## Coding Style

- 4 spaces for indentation rather than tabs
- 80 character line length up to 120 when needed, but try to keep it down to 80
  chars per line
- Use JuliaFormatter.jl to format your files
