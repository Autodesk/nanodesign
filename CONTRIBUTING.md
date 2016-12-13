# Contributing to the Nanodesign Python Package

# Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Project Roles](#project-roles)
    - [Maintainers](#maintainers)
    - [Contributors](#contributors)
- [Important Details](#important-details)
    - [Branching Model](#branching-model)
    - [Coding Style Guide](#coding-style-guide)
    - [Tests](#tests)
    - [Documentation](#documentation)
    - [Feature Requests and Questions](#feature-requests-and-questions)
- [What can I contribute?](#what-can-i-contribute)
    - [General Contributions](#general-contributions)
    - [Important note about pull requests](#important-note-about-pull-requests)
    - [Design and cleanup proposals](#design-and-cleanup-proposals)
- [Submitting an Issue or Pull Request](#submitting-an-issue-or-pull-request)
    - [Timing](#timing)
    - [Submitting an Issue](#submitting-an-issue)
    - [Submitting Your Pull Request](#submitting-your-pull-request)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->
<!-- to generate: npm install doctoc: doctoc --gitlab --maxlevel 3 CONTRIBUTING.md-->


## Project Roles
### Maintainers

Maintainers are responsible for responding to pull requests and issues, as well as guiding the direction of the project. 

We currently have the following maintainers:

* Joseph Schaeffer <joseph.schaeffer@autodesk.com> - Project Lead, Maintainer
* Dave Parker <dave.parker@autodesk.com> - Maintainer

If you've established yourself as an impactful contributor for the project, and are willing take on the extra work, we'd love to have your help maintaining it! Email the maintainers list at `nanodesign_maintainers@autodesk.com` for details.

### Contributors

We encourage everyone to contribute to this project! You can become a contributor by forking this project and submitting pull requests. You can also contribute by reporting bugs in the Issue tracker, helping review pull requests, participate in discussions on the forums about the package, and more. 

This document exists to help guide people in what they can contribute and how they can contribute.

----

## Important Details
### Branching Model

We use the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model. This means that all pull requests should be onto the *dev* branch. The *master* branch should always have a version tag and be the latest released version. Upcoming releases will be found in branches with the *release* prefix. 

We prefer all pull requests come from forks of the repository, rather than having maintainers and contributors create branches within the main repository. 

### Coding Style Guide

We intend to follow the [PEP-8](https://www.python.org/dev/peps/pep-0008/) style guide. The [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html) is also a good reference, but we don't follow it as closely. If you find sections of code that don't meet the style recommendations in spirit as well as letter, please feel free to contribute in that area via raising an issue and then possibly filing a pull request addressing the issue.

### Tests

We currently have a variety of end-to-end tests which test many of the scripts included with the package. New contributions should ensure that these end-to-end tests are still functioning correctly with any code changes. These tests can be found within the [`tests/`](tests/) subdirectory, and there are specific subdirectories there for each of the main script files, currently `converter`, `visualizer`, and `stapler`. 

We are also intending to have unit tests for individual modules. These will reside in the same `tests/` directory structure rather than directly with the code. These should use the pytest module for collecting and running the tests. Contributions and suggestions are very welcome for how we should best implement a unit testing strategy.

Once we have solidified our testing pieces, we will be enforcing the requirement that *all* pull requests must include associated tests for the components affected. These tests should eventually grow into something which we can automate via Continuous Integration tools to streamline the building and verification process for releases.

### Documentation

All modules, classes, functions, etc, should be appropriately documented. In addition, we are planning on having both automated API documentation as well as curated documentation available. The [`docs/`](docs/) subdirectory contains the basic code to generate this documentation, as well as the start of curated documentation guides. This area is still a work in progress, and contributions are always welcome for improving the documentation.

For a walkthrough of using some of the existing features, check out the [`scripts/`](scripts/) directory for three examples of using the package to perform different operations:
* [`converter.py`](scripts/converter.py): Convert between different input and output formats.
* [`stapler.py`](scripts/stapler.py): Given a Cadnano input file, automatically create staples given the scaffold structure. Has options for removing existing staples, removing a subset of staples, and using existing staples and only breaking those.
* [`vis.py`](scripts/vis.py): Visualize a structure, uses PyOpenGL for the 3D scene.

### Feature Requests and Questions

If you have feature requests or questions, please use our [**Forums**](https://forum.bionano.autodesk.com/c/Nano-Design/nanodesign-python).

Alternately, you can file an issue on Github, but the forums are preferred. If you have an opinion on this guideline, let us know!

-------

## What can I contribute?

If you're interested in getting started, here are some general contributions that are always welcome. Once we've launched, we'll start maintaining a wishlist of specific ideas on the [forums](https://forum.bionano.autodesk.com/c/Nano-Design/nanodesign-python) and possibly in the github wiki for the project.

### General Contributions

**Tests**: This is one of the easiest ways to get started - we currently have a lot of end-to-end tests, but more specific end-to-end tests and new unit tests are very welcome!

 * **Contribute unit tests** - test that Nanodesign's functions, classes and methods work correctly in a variety of situations, and that they report errors when appropriate. We're planning on using `pytest` for this, and if you have feedback on how to best implement unit tests for this package, we'd love to hear it!
 * **Contribute interesting sample files** - We can always use sample files that test the boundaries of the input parser(s). You'll see several such files already in [`tests/samples/`](tests/samples/), and if you find a file that causes issues, adding it to the sample files so we can run tests on it is appreciated.
 * **Contribute end-to-end tests** - We currently have end-to-end tests that use our basic scripts. If you have more complicated tests that can run against specific workflows and provide useful error handling, those could be a good contribution to the repository!

**Examples**: We would love to have more example scripts showing off the functionality. See the [`scripts/`](scripts/) directory for some existing examples.

 * **Specific tasks** - put together a script that runs a specific task on an input. For example, something that loads a cadnano file and finds all domains that have a particular sequence. 
 * **Interesting use cases** - a script that runs some complicated construction or algorithm, such as building an arbitrary length N-helix bundle.

**Bug fixes:** Found a typo in the code? Found that a function fails under certain conditions? Know how to fix it? Great! Go for it. Please do [open an issue](https://github.com/autodesk/nanodesign/issues) so that we know you're working on it, and submit a pull request when you're ready.

**Features:** We'd love to have new features be added. Please [open an issue](https://github.com/autodesk/nanodesign/issues) so that we know what you're working on.

### Important note about pull requests

All PRs should be documented as [GitHub issues](https://github.com/autodesk/nanodesign/issues), ideally BEFORE you start working on them.

### Design and cleanup proposals

We understand this project is a work in progress, and an important part of that is improving the API, cleaning up the code, adding tests, and writing more documentation. We would love for help on these, even if you don't do any programming! 

For example:
 * You could describe how a user should be able to construct a new DNA structure from scratch. Do they describe it via contiguous helices? Do they give the topology only? How do we construct a 3D model given these constraints?
 * You can also propose redesigns or package changes - for example you might want to move the stapler algorithm directly into the core package, and think it should be `nanodesign.algorithms.stapler`, and provide access to specific functions rather than the general `Stapler` class. 

To get started, as always: [open an issue](https://github.com/autodesk/nanodesign/issues). Let us know what you think, and we'll figure out a way to add those into the wishlist or design documents!


------

## Submitting an Issue or Pull Request

### Timing

We will attempt to address all issues and pull requests within one week. It may a bit longer before pull requests are actually merged, as they must be inspected and tested. 

### Submitting an Issue

If the Nanodesign python package isn't working like you expect, please open a new issue! We appreciate any effort you can make to avoid reporting duplicate issues, but please err on the side of reporting the bug if you're not sure.

Providing the following information will increase the chances of your issue being dealt with quickly:

* **Overview of the Issue** - Please describe the issue, and include any relevant exception messages or screenshots.
* **Environment** - Include the relevant output of `pip freeze` as well as your system and python version info.
* **Help us reproduce the issue** - Please include code that will help us reproduce the issue. For complex situations, attach a python script exhibiting the behavior.
* **Related Issues** - Please link to other issues in this project (or even other projects) that appear to be related 

### Submitting Your Pull Request

Before you submit your pull request consider the following guidelines:

* Search GitHub for an open or closed Pull Request or Issue that relates to your submission. You don't want to duplicate effort.
* Make your changes in a new git branch:

     ```shell
     git checkout -b my-pr-branch [working-branch-name]
     ```

* Create your patch.
* Commit your changes using a descriptive commit message.

     ```shell
     git add <changed files>
     git commit
     ```
  Note: we do not recommend using the `-a` command line option here, as while it will automatically "add" and "rm" edited files, it can frequently add things you do not wish to have included, such as testing files not in the `.gitignore`. Specifying the exact files you want included is preferable.
  
  We recommend (but do not require) a descriptive commit/PR message that includes the following three things:
  - Short descriptive title (~50 chars), answering what the commit/PR does.
  - Why is a change needed?  Explain why this fix was needed. This can be longer and should give the motivation for this change. Should reference issues or discussions as needed.
  - How did we implement this change? Explain how it works and why it answers the need.

  If you're interested in automating this type of message template, or for more background, read [this link](https://robots.thoughtbot.com/better-commit-messages-with-a-gitmessage-template). 

* Push your branch to GitHub:

    ```shell
    git push origin my-fix-branch
    ```

* In GitHub, send a pull request to `nanodesign:dev`
* Before any request is merged, you'll need to agree to the contribution boilerplate. This should occur automatically when you file the pull request, but you can email us at `nanodesign_maintainers@autodesk.com` for more detail if needed. 
