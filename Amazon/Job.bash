#!/bin/bash 

#clean temp file for memory utility
rm -rf /var/tmp/aws-mon.bak
mv /var/tmp/aws-mon /var/tmp/aws-mon.bak

cd /dev

#clone
git clone --depth 1  git@github.com:bw4sz/WhaleShape.git

cd WhaleShape||sudo halt

#make new branch
#name it the instance ID
iid=$(ec2metadata --instance-id)

git checkout -b $iid

#render script
Rscript -e "rmarkdown::render('DynamicForaging.Rmd')" &> run.txt

#push results
git add --all
git commit -m "ec2 run complete"
git push -u origin $iid

#kill instance
sudo halt
