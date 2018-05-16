cat workshop_trajectories.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd
echo "rmarkdown::render('answers.Rmd', clean=FALSE)" | R --slave --no-save
mv answers.knit.md answers.md
