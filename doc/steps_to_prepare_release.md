# Steps to prepare a release of Salmon
-----

 1. Tag corresponding commit of RapMap so that it can be stably pulled in for source builds.
 2. Alter `fetchRapMap.sh` to fetch the corresponding tagged version (and update the sha256 sum).
 3. Bump salmon version in `include/SalmonConfig.hpp`, then rebuild and run the `bump_version.sh` script.
 4. Ensure that everything builds cleanly on Linux (taken care of by CI) and OSX.
 5. Merge the develop branch changes into master.
 6. Tag the Salmon release with a new version number.
 7. Update the docker tag and build an image for docker hub.
 8. Bump the Bioconda version and build a new Bioconda release.
 9. Add release notes for the tagged master version.
 10. Upload the pre-compiled linux binary (from the CI server) to GitHub.
 11. (not technically part of release) Reset the relevant changes (steps 1,2) on the develop branch so they now point to a non-tagged RapMap.
