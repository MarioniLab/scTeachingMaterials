cat workshop_trajectories.Rmd | grep -v "^#A#" | sed "s/^#H# //" > practical.Rmd
cat workshop_trajectories.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd
rm workshop_trajectories.Rmd

rm answers.md
rm test.sh
