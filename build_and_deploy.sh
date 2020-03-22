#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Then we build and deploy the sphinx / doxygen documentation
SOURCE_BRANCH="development"

# Pull requests and commits to other branches shouldn't try to deploy
if [ "$TRAVIS_PULL_REQUEST" != "false" -o "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" ]; then
    echo "Skipping deploy."
    exit 0
fi

# Save some useful information
REPO=`git config remote.origin.url`
#SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`

# sphinx
#
# run breathe and clean up
cd Docs/sphinx_documentation
#
#breathe-apidoc --o source ../../out/docs_xml/doxygen/ -g class,file
#python make_api.py
#
echo "Build the Sphinx documentation for incflo."
make SPHINX_BUILD="python -msphinx" latexpdf
mv build/latex/incflo.pdf source/ 
make SPHINX_BUILD="python -msphinx" html &> make_source_html.out
#
# Start ssh-agent
#openssl aes-256-cbc -K $encrypted_11fd376b52bf_key -iv $encrypted_11fd376b52bf_iv -in deploy_rsa.enc -out deploy_rsa -d
#chmod 600 ./deploy_rsa
#eval `ssh-agent -s`
#ssh-add ./deploy_rsa

# clone document
DOC_SSH_REPO="git@github.com:AMReX-Codes/AMReX-Codes.github.io.git"
git clone $DOC_SSH_REPO AMReX-Codes.github.io
cd AMReX-Codes.github.io/incflo
#
# clean out existing contents
git rm -r docs_html tutorials_html docs_xml
mkdir     docs_html tutorials_html docs_xml
#
# add sphinx
cp -rp ../../Docs/sphinx_documentation/build/html/* docs_html/
cp -rp ../../Docs/sphinx_tutorials/build/html/* tutorials_html/
#
# Now let's go have some fun with the cloned repo
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

if git diff-index --quiet HEAD; then
    exit 0
fi

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add docs_html tutorials_html docs_xml
git commit -m "Deploy to GitHub Pages: ${SHA}" || true

git push $DOC_SSH_REPO || true
ssh-agent -k
cd ..
