# ubiDash

Interactive dashboard to query multi-omic datasets from human tissues and/or cell lines, with special focus on E3 ubiquitin ligases of the Cullin RING Ligase (CRL) family.

## Internal Git teamwork protocol

Before starting your work, make sure to setup the remote repo appropriately. Then you can follow step referenced [here](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).

1.  Before creating a new branch, pull the changes from upstream. Your master needs to be up to date.

    -   `$ git pull` or `$ git pull --rebase`

2.  Create a branch on your local machine and switch into this branch.

    -   `$ git checkout -b [name_of_your_new_branch]` or
    -   just create a branch by `$ git branch [name_of_your_new_branch]`

3.  Start working on your local branch. After satisfied with changes made, follow the workflow below:

    -   `$ git add .` to stage all changes

    -   `$ git commit -m "add descriptive message here"` to commit all changes

    -   `$ git push origin [name_of_your_new_branch]` to push all committed changes to the remote branch

    -   You can see all the branches created by using: `$ git branch -a`

4.  To merge your branch to the master branch:

    -   To update original branch when branch from official repository has been updated: `$ git fetch [name_of_your_remote]`
    -   Merge changes: `$ git merge [name_of_your_remote]/[name_of_your_new_branch]`
    -   Delete a branch on your local filesystem: `$ git branch -d [name_of_your_new_branch]`
    -   Delete the branch on github: `$ git push origin :[name_of_your_new_branch]`
