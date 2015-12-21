#!/bin/sh
cd ./src/
SESSION="netpop"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n net34
tmux send-keys -t $SESSION:1 'vim utils.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe utils.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe network.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe network.c++' C-m

tmux split-window -h
tmux send-keys -t $SESSION:1 'vim net3.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net3.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net4.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net4.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net4-homog.c++' C-m
tmux select-pane -t 0

tmux new-window -t $SESSION:2 -k -n dendlatt
tmux send-keys -t $SESSION:2 'vim randnet_dendritic.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe randnet_dendritic.c++' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe randnet_lattice.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe randnet_lattice.c++' C-m

tmux split-window -h
tmux send-keys -t $SESSION:2 'vim utils_dendlatt.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe utils_dendlatt.c++' C-m

cd ../
tmux new-window -t $SESSION:3 -n makefile
tmux select-window -t $SESSION:3
tmux send-keys -t $SESSION:3 'vim .travis.yml' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe tmux-start.bash' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe README.md' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe CMakeLists.txt' C-m

tmux split-window -h
tmux select-pane -t 0

cd ./src/
tmux new-window -t $SESSION:4 -n R
tmux select-window -t $SESSION:4
tmux send-keys -t $SESSION:4 'vim net2.r' C-m
tmux send-keys -t $SESSION:4 ':' 'tabe aaajunk.r' C-m

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
