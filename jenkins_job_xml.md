# Retrieve and/or backup a Jenkins job

It may be possible to run the Jenkins backup manager as a non-admin, but
retrieving the xml file for a job requires access to the institutional W-13 user
account ``toolbox-jenkins``. Access is reserved for the W-13 DevOps team.

To get a copy of a Jenkins job xml file for revision control or backup, contact
a member of the W-13 DevOps team. What follows is the manual used by the DevOps
team.

1. From the [Jenkins landing/main page](https://toolbox-jenkins.lanl.gov/)
   select "Manage Jenkins" in the left sidebar.

2. Scroll all the way down to select "Backup Manager" in the "Uncategorized"
   section.

3. You can double check the settings

   * Backup configuration

     * Hudson root directory ``/local/jenkins-homes/toolbox-jenkins``
     * Backup directory: ``/home/toolbox-jenkins/backups``
     * Format: ``tar.bz2``
     * File name template: ``bakup_@date@.@extension@``
     * Custom exclusions: ``archive``
     * Verbose mode: checked
     * Configuration files (.xml) only: unchecked
     * No shutdown: checked

   * Backup content

     * Backup job workspace: unchecked **NOTE: checking this box will vastly
       increase the backup file's disk space requirements! Use with caution.**
     * Backup builds history: checked **Note: these are relatively small build
       meta data files about success/failure.**
     * Backup maven artifacts archives: checked **Note: these file sizes depend
       entirely on what a project owner decides to archive to ninetails after
       the build. This option could be large if someone decided to copy anything
       other than pytest results xml files or small text files. If everyone is
       following good practices, these should be small.**
     * Backup fingerprints: checked

4. Press the "Backup Hudson Configuration" and wait

   * If you uncheck "backup configuration file only" the backup can take several
     hours because it will save build, configuration, and artifacts for
     EVERY build of EVERY job. With that option unchecked, the backup often
     still requires 40 minutes.
   * If you click away from the "Backup manager log" webpage, you won't be able
     to get back to it. You'll be waiting blind for the job to finish. You can
     check progress by watching the backup file build and change size at

```
toolbox-jenkins@ninetails:~/backups$ pwd
/home/toolbox-jenkins/backups
toolbox-jenkins@ninetails:~/backups$ ll -h
total 1.3G
-rw-r--r--. 1 toolbox-jenkins toolbox-jenkins 278M Apr  7  2020
bakup_20200407_0939.tar.bz2
-rw-r--r--. 1 toolbox-jenkins toolbox-jenkins 284M Apr 30 13:37
bakup_20200430_1312.tar.bz2
-rw-r--r--. 1 toolbox-jenkins toolbox-jenkins 284M Jun  8 14:16
bakup_20200608_1352.tar.bz2
-rw-r--r--. 1 toolbox-jenkins toolbox-jenkins 289M Jun 22 15:20
bakup_20200622_1453.tar.bz2
-rw-r--r--. 1 toolbox-jenkins toolbox-jenkins 183M Oct 22 11:32
bakup_20201022_1119.tar.bz2
```
with a command like
```
toolbox-jenkins@ninetails:~/backups$ watch ls -lh bakup_20201022_1119.tar.bz2
```
When it stops changing size for several minutes, the backup is probably
complete.
