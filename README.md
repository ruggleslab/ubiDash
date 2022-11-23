# ubiDash

Interactive dashboard to query multi-omic datasets from human tissues and/or cell lines, with special focus on E3 ubiquitin ligases of the Cullin RING Ligase (CRL) family.

## Internal Git teamwork protocol

Before starting your work, make sure to setup the remote repo appropriately. Then you can follow step referenced [here](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).

1.  Before creating a new branch, pull the changes from upstream. Your master needs to be up to date.

    -   `$ git pull` or `$ git pull --rebase`

2.  Create a branch on your local machine and switch into this branch. (Refs: [1](https://www.git-tower.com/learn/git/commands/git-checkout), [2](https://www.nobledesktop.com/learn/git/git-branches))

    -   To create and switch to new branch use `$ git checkout -b [name_of_your_new_branch]` or
    -   just create a branch by `$ git branch [name_of_your_new_branch]`
    -   To set remote branch tracking on new branch `$ git checkout -b <new-branch> --track <remote-branch>`
    -   To list all existing branches use `$ git branch -a`, `$ git branch` for local and `$ git branch -r` for remote repos.
    -   To display more branch details use `$ git branch -vv` or `$ git branch -vva`

3.  Start working on your local branch. After satisfied with changes made, follow the workflow below:

    -   `$ git add .` to stage all changes

    -   `$ git commit -m "add descriptive message here"` to commit all changes

    -   `$ git push origin [name_of_your_new_branch]` to push all committed changes to the remote branch

    -   You can see all the branches created by using: `$ git branch -a`

4.  Merging branches to master

    -   Make sure your working tree is clean (`$ git status`)

    -   First, must checkout branch you want to merge another branch into (changes will be merged into this branch). If not currently on the desired branch run `$ git checkout master`

        -   **NOTE:** Replace `master` with another branch name as needed.

    -   Now, can merge another branch into current branch by `$ git merge my-branch-name`

        -   **NOTE:** When you merge, there may be a conflict. Refer to [**Handling Merge Conflicts**](https://www.nobledesktop.com/learn/git/merge-conflicts) to learn what to do.
        -   You can also look into this [link](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-comparing-branches-in-pull-requests) to learn how to compare between branches and branch versions.

5.  Deleting existing branches you no longer need

    -   To delete a remote branch use `$ git push origin --delete my-brnach-name`

    -   To delete a local branch use `$ git branch -d my-branch-name` or `$ git branch -D my-branch-name`

        -   **NOTE:** The -d option only deletes the branch if it has already been merged. The -D option is a shortcut for \--delete \--force, which deletes the branch irrespective of its merged status.
