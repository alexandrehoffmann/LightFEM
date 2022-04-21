#!/bin/bash

for file in u_*.dat; do
   gnuplot -e "set terminal jpeg large font arial size 1024,768; set view 0,0; set palette defined (-0.05 \"blue\", 0.0 \"white\", 0.05 \"red\"); set cbrange[-0.05:0.05]; splot '$file' u 1:2:3 w pm3d" > pic_$(echo "$file" | cut -f 1 -d '.').jpeg
done

ffmpeg -s 1024x768 -i pic_u_%d.jpeg -crf 15 u.mpeg
rm *.jpeg

for file in sem_u_*.dat; do
   gnuplot -e "set terminal jpeg large font arial size 1024,768; set view 0,0; set palette defined (-0.05 \"blue\", 0.0 \"white\", 0.05 \"red\"); set cbrange[-0.05:0.05]; splot '$file' u 1:2:3 w pm3d" > pic_$(echo "$file" | cut -f 1 -d '.').jpeg
done

ffmpeg -s 1024x768 -i pic_sem_u_%d.jpeg -crf 15 sem_u.mpeg
rm *.jpeg

ffmpeg -i sem_u.mpeg -i u.mpeg \
 -filter_complex \
    "[0:v]pad=iw*2:ih[int]; \
     [int][1:v]overlay=W/2:0[vid]" \
-map "[vid]" \
-c:v libx264 -crf 15 \
sem_vs_fem.mp4
