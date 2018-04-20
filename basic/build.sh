cat workshop_kolod.Rmd | grep -v "^#A#" | sed "s/^#H# //" > practical.Rmd
cat workshop_kolod.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd

