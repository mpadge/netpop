#!/bin/sh
cd ./src/
SESSION="netpop"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n main
tmux send-keys -t $SESSION:1 'vim utils.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe utils.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net3.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net3.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net4.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe net4.c++' C-m

cd ../
tmux split-window -h
tmux select-pane -t 0

cd ./src/
tmux new-window -t $SESSION:2 -n R
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim net2.r' C-m

cd ../
tmux new-window -t $SESSION:3 -n makefile
tmux select-window -t $SESSION:3
tmux send-keys -t $SESSION:3 'vim makefile' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe tmux-start.bash' C-m

tmux split-window -h
tmux select-pane -t 0

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
