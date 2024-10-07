ffmpeg -i GFP3.avi -i vis-animated-smaller.mp4 \
-filter_complex "[1:v]setpts=(8/33.6)*PTS[v2]; \
[0:v]hflip,scale=640:360[v1]; [v2]scale=640:360[v2_resized]; \
[v1][v2_resized]hstack=inputs=2;" \
     -c:v libx264 output.mp4

ffmpeg -i output.mp4 -filter:v "setpts=PTS*3" -c:v libx264 output_slower.mp4
