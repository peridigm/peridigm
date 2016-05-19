# Contributing to Peridigm

## How To Open a Peridigm Pull Request (and get it merged)

The following commands are used by Peridigm contributors to set up a fork of Peridigm's Git repository on GitHub, create a new branch, and create a GitHub pull request ("PR") of the changes in that branch.

## Set up your own fork of the Peridigm repository

 1. [Fork the peridigm/peridigm repository on GitHub](https://github.com/peridigm/peridigm/fork).
      * This creates a personal remote repository that you can push to. This is needed because only Peridigm maintainers have push access to the main repositories.
 1. Clone your forked repository with `git clone https://github.com/<YOUR_GIHUB_USER_NAME>/peridigm.git`
 1. Change to the directory containing your Peridigm installation with `cd peridigm`
 1. Add the main Peridigm repository as an upstream remote with `git remote add upstream https://github.com/peridigm/peridigm.git`.

## Create your pull request from a new branch

To make a new branch and submit it for review, create a GitHub pull request with the following steps:

 1. Check out the `master` branch with `git checkout master`.
 2. Retrieve new changes to the `master` branch from the main Peridigm repository with `git fetch upstream master` followed by `git rebase upstream/master`.
 3. Create a new branch from the latest `master` branch with `git checkout -b <YOUR_BRANCH_NAME> origin/master`.
 4. Make your changes, bug fixes, or features additions to the Peridigm code.
 5. Ensure that your new changes do not affect any existing code that has tests.  [Build Peridigm](https://github.com/peridigm/peridigm/blob/master/doc/BuildingPeridigm.md) with your new changes.  After successfully compiling run `make test` and verify that all tests pass.
    * If your pull request is anything other than a bug-fix, it's unlikely it will be accepted into the Peridigm repository without sufficient tests of the new   feature, so please consider writing your own unit, regression, and/or verification tests and include them in the PR
 6. Make a separate commit for each meaningful change with `git add` and `git commit`.
 7. Perform a squash merge from you new branch back to the master branch
    1. `git checkout master`
    2. `git merge --squash <YOUR_BRANCH_NAME>`
    3. `git commit`
 8. Upload your new commits to the branch on your fork with `git push origin master`.
 9. Go to https://github.com/peridigm/peridigm and create a pull request to request review and merging of the commits in your pushed branch. Explain why the change is needed and, if fixing a bug, how to reproduce the bug. 
10. Await feedback or a merge from Peridigm's maintainers. We typically respond to all PRs within a couple days, but it may take up to a week, depending on the maintainers' workload.
11. Thank you!

## Following up

To respond well to feedback:

 1. Ask for clarification of anything you don't understand and for help with anything you don't know how to do.
 2. Post a comment on your pull request if you've provided all the requested changes/information and it hasn't been merged after a week. Post a comment on your pull request if you're stuck and need help.
    * A `needs response` label on a PR means that the Peridigm maintainers need you to respond to previous comments.
 3. Keep discussion in the pull request unless requested otherwise (i.e. do not email maintainers privately).
 4. Do not continue discussion in closed pull requests.
 5. Do not argue with Peridigm maintainers. You may disagree but unless they change their mind, please implement what they request. Ultimately they control what is included in Peridigm, as they have to support any changes that are made.

To make changes based on feedback:

 1. Check out your branch again with `git checkout <YOUR_BRANCH_NAME>`.
 2. Make any requested changes and commit them with `git add` and `git commit`.
 3. Squash new commits into one commit per formula with `git rebase --interactive origin/master`.
 4. Push to your remote fork's branch and the pull request with `git push --force`.

If you are working on a PR for a small change that was originally just one commit, `git commit --amend` is a convenient way of keeping your commits squashed as you go.

Once all feedback has been addressed and if it's a change we want to include, then we'll add your commit to Peridigm. Note that the PR status may show up as "Closed" instead of "Merged" because of the way we merge contributions. Don't worry: you will still get author credit in the actual merged commit.

Well done, you are now a Peridigm contributor!

---

**IMPORTANT:** By submitting a patch, you agree to allow the project owners to license your work under the terms of the [Peridigm License](https://peridigm.sandia.gov/content/license).
